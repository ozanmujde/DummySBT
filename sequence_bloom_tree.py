import numpy as np
from rawhash import RawHash
from collections import defaultdict

class SequenceBloomTree:
    """
    N-ary Sequence Bloom Tree implementation for efficient sequence comparison
    """
    def __init__(self, k_mer_size=20, bloom_size=1048576, hash_count=3, branching_factor=4):
        """
        Initialize Sequence Bloom Tree
        
        :param k_mer_size: Size of k-mers (default 20 as per HowDeSBT)
        :param bloom_size: Size of Bloom filter (must match RawHash)
        :param hash_count: Number of hash functions (must match RawHash)
        :param branching_factor: Number of children per internal node
        """
        self.k_mer_size = k_mer_size
        self.bloom_size = bloom_size
        self.hash_count = hash_count
        self.branching_factor = branching_factor
        self.tree = None
        self.leaf_map = {}  # Map leaf index to sequence index
        
    def _initialize_tree(self, num_sequences):
        """
        Create an n-ary Bloom filter tree
        
        :param num_sequences: Number of sequences
        :return: Dictionary containing tree structure
        """
        # Calculate number of internal nodes needed for n-ary tree
        num_leaves = num_sequences
        height = max(1, int(np.ceil(np.log(num_leaves) / np.log(self.branching_factor))))
        num_internal = int((self.branching_factor ** height - 1) / (self.branching_factor - 1))
        total_nodes = num_internal + num_leaves
        
        # Initialize nodes with zeros
        nodes = [np.zeros(self.bloom_size, dtype=np.bool_) for _ in range(total_nodes)]
        
        return {
            'nodes': nodes,
            'num_internal': num_internal,
            'total_nodes': total_nodes,
            'height': height
        }
        
    def _get_children_idx(self, node_idx):
        """Get children indices for a given node"""
        start = node_idx * self.branching_factor + 1
        return range(start, min(start + self.branching_factor, len(self.tree['nodes'])))
        
    def _get_parent_idx(self, node_idx):
        """Get parent index for a given node"""
        if node_idx <= 0:
            return -1
        return (node_idx - 1) // self.branching_factor
        
    def insert_signatures(self, signatures, sequences=None):
        """Insert signatures into the tree"""
        if not signatures:
            return
            
        # Initialize tree if needed
        if not self.tree:
            self.tree = self._initialize_tree(len(signatures))
            
        # Clear existing mappings
        self.leaf_map.clear()
        
        # Insert signatures as leaf nodes
        for i, sig in enumerate(signatures):
            leaf_idx = self.tree['num_internal'] + i
            if leaf_idx < len(self.tree['nodes']):
                self.tree['nodes'][leaf_idx] = sig.copy()
                self.leaf_map[leaf_idx] = i
                
                # Update parent nodes using simple OR
                current = self._get_parent_idx(leaf_idx)
                while current >= 0:
                    # Get all children signatures
                    children = [self.tree['nodes'][idx] for idx in self._get_children_idx(current)
                              if idx < len(self.tree['nodes'])]
                    if children:
                        # Combine using OR
                        combined = children[0].copy()
                        for child in children[1:]:
                            combined |= child
                        self.tree['nodes'][current] = combined
                    current = self._get_parent_idx(current)
                    
    def query(self, query_signature, threshold=0.7, expected_matches=None):
        """Query the tree for matching sequences"""
        if not self.tree:
            raise ValueError("Tree not initialized")
            
        def containment_score(query, node):
            """Calculate containment score"""
            intersection = np.sum(query & node)
            query_ones = np.sum(query)
            return intersection / query_ones if query_ones > 0 else 0
            
        def check_node(node_idx, query_sig):
            """Check node and its children for matches"""
            if node_idx >= len(self.tree['nodes']):
                return set()
                
            node = self.tree['nodes'][node_idx]
            score = containment_score(query_sig, node)
            
            # Early exit if score is too low
            if score < threshold * 0.5:  # More lenient for internal nodes
                return set()
                
            # If leaf node, check final threshold
            if node_idx >= self.tree['num_internal']:
                return {node_idx} if score >= threshold else set()
                
            # Check all children
            matches = set()
            for child_idx in self._get_children_idx(node_idx):
                matches.update(check_node(child_idx, query_sig))
            return matches
            
        # Find matching leaves
        matching_leaves = check_node(0, query_signature)
        
        # Convert leaf indices to sequence indices
        matched_sequences = {self.leaf_map[leaf_idx] for leaf_idx in matching_leaves
                           if leaf_idx in self.leaf_map}
        
        # Calculate metrics
        if expected_matches is not None:
            expected_set = set(expected_matches)
            true_positives = len(matched_sequences & expected_set)
            false_positives = len(matched_sequences - expected_set)
            false_negatives = len(expected_set - matched_sequences)
        else:
            # If no expected matches provided, assume query signature should match itself
            true_positives = len(matched_sequences)
            false_positives = 0
            false_negatives = 0 if matched_sequences else 1
        
        # Calculate final metrics
        precision = true_positives / (true_positives + false_positives) if true_positives + false_positives > 0 else 0
        recall = true_positives / (true_positives + false_negatives) if true_positives + false_negatives > 0 else 0
        f1_score = 2 * (precision * recall) / (precision + recall) if precision + recall > 0 else 0
        
        return {
            'matches': list(matched_sequences),
            'metrics': {
                'true_positives': true_positives,
                'false_positives': false_positives,
                'false_negatives': false_negatives,
                'precision': precision,
                'recall': recall,
                'f1_score': f1_score
            }
        }

# Example usage
if __name__ == "__main__":
    # Create dummy signatures for demonstration
    sbt = SequenceBloomTree()
    dummy_sequences = ['ACGT' * 25 for _ in range(10)]  # 100bp sequences
    dummy_signatures = [
        np.random.randint(0, 2, size=1 << 20, dtype=np.bool_) 
        for _ in range(10)
    ]
    
    sbt.insert_signatures(dummy_signatures, dummy_sequences)
    
    # Query with a random signature
    query_sig = np.random.randint(0, 2, size=1 << 20, dtype=np.bool_)
    results = sbt.query(query_sig)
    
    print(f"Query results: {results}")
