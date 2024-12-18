import json
import mmh3
import numpy as np
import random
from typing import List, Set, Dict, Tuple

class RawHash:
    """
    Custom hash method for processing sequence data with improved Bloom filter generation
    """
    def __init__(self, k_mer_size=20, bloom_size=1048576, hash_count=3):
        """
        Initialize RawHash
        
        :param k_mer_size: Size of k-mers
        :param bloom_size: Size of Bloom filter (fixed to avoid inconsistencies)
        :param hash_count: Number of hash functions
        """
        self.k_mer_size = k_mer_size
        self.bloom_size = bloom_size
        self.initial_hash_count = hash_count  # Store initial hash count
        self.hash_count = hash_count
        self.kmer_cache = {}  # Cache for k-mer hashes
        self.max_density = 0.4  # Maximum allowed density
        
    def generate_kmers(self, sequence: str) -> List[str]:
        """
        Generate k-mers from a given sequence with error checking
        
        :param sequence: Input sequence
        :return: List of valid k-mers
        """
        if not sequence or len(sequence) < self.k_mer_size:
            return []
            
        # Convert to uppercase and validate
        sequence = sequence.upper()
        valid_bases = set('ACGT')
        
        # Pre-validate the entire sequence
        if not all(base in valid_bases for base in sequence):
            # Clean the sequence first
            sequence = ''.join(base for base in sequence if base in valid_bases)
            if len(sequence) < self.k_mer_size:
                return []
        
        # Use set for unique k-mers
        kmers = set()
        for i in range(len(sequence) - self.k_mer_size + 1):
            kmer = sequence[i:i + self.k_mer_size]
            kmers.add(kmer)
        
        return list(kmers)
        
    def hash_kmer(self, kmer: str, seed: int = 0) -> int:
        """
        Hash a k-mer using cached MurmurHash3 with a seed
        
        :param kmer: K-mer string to hash
        :param seed: Seed for hashing
        :return: Hash value
        """
        cache_key = (kmer, seed)
        if cache_key not in self.kmer_cache:
            # Convert k-mer to bytes and include position information
            kmer_bytes = kmer.encode('utf-8')
            # Use both 32-bit and 64-bit hashes for better distribution
            h1 = mmh3.hash64(kmer_bytes, seed=seed)[0]
            h2 = mmh3.hash(kmer_bytes, seed=seed)  # 32-bit hash
            # Combine hashes
            self.kmer_cache[cache_key] = h1 ^ ((h2 & 0xFFFFFFFF) << 32)
        return self.kmer_cache[cache_key]
        
    def calculate_hash_positions(self, kmer: str) -> List[int]:
        """
        Calculate all hash positions for a k-mer
        
        :param kmer: Input k-mer
        :return: List of hash positions
        """
        h1 = abs(self.hash_kmer(kmer, 0))
        h2 = abs(self.hash_kmer(kmer, len(kmer)))
        h3 = abs(self.hash_kmer(kmer, h1 % 100))
        
        positions = []
        for i in range(self.hash_count):
            pos = (h1 + (i * h2) + ((i * i) * h3)) % self.bloom_size
            positions.append(pos)
            
        return positions
        
    def adjust_for_density(self, kmers: List[str], hash_positions: Dict[str, List[int]]) -> Tuple[List[str], int]:
        """
        Adjust k-mers or hash count to meet density requirements
        
        :param kmers: List of k-mers
        :param hash_positions: Dictionary of k-mer to hash positions
        :return: Tuple of (adjusted k-mers, adjusted hash count)
        """
        # Count unique positions
        all_positions = set()
        for positions in hash_positions.values():
            all_positions.update(positions)
            
        density = len(all_positions) / self.bloom_size
        
        if density <= self.max_density:
            return kmers, self.hash_count
            
        # Try reducing hash count first
        min_hash_count = max(1, int(self.max_density * self.bloom_size / len(kmers)))
        if min_hash_count < self.hash_count:
            return kmers, min_hash_count
            
        # If still too dense, sample k-mers
        target_kmers = int(self.max_density * self.bloom_size / self.hash_count)
        sampled_kmers = random.sample(kmers, min(target_kmers, len(kmers)))
        return sampled_kmers, self.hash_count
        
    def process_read(self, read: Dict[str, str]) -> np.ndarray:
        """
        Process a single read and generate its hash signature with improved density control
        
        :param read: Dictionary containing read data
        :return: Bloom filter signature for the read
        """
        sequence = read['sequence']
        
        # Generate k-mers
        kmers = self.generate_kmers(sequence)
        if not kmers:
            return np.zeros(self.bloom_size, dtype=np.bool_)
        
        # Reset hash count to initial value
        self.hash_count = self.initial_hash_count
        
        # First calculate all hash positions
        hash_positions = {kmer: self.calculate_hash_positions(kmer) for kmer in kmers}
        
        # Adjust for density
        kmers, self.hash_count = self.adjust_for_density(kmers, hash_positions)
        
        # Create and populate Bloom filter
        bloom_filter = np.zeros(self.bloom_size, dtype=np.bool_)
        for kmer in kmers:
            positions = hash_positions[kmer][:self.hash_count]  # Use only needed positions
            for pos in positions:
                bloom_filter[pos] = True
                
        return bloom_filter
        
    def process_file(self, filepath: str) -> List[np.ndarray]:
        """
        Process a JSON file containing reads
        
        :param filepath: Path to JSON file
        :return: List of Bloom filter signatures
        """
        with open(filepath, 'r') as f:
            data = json.load(f)
            return [self.process_read(read) for read in data]

# Example usage
if __name__ == "__main__":
    rawhash = RawHash()
    print("RawHash initialized with improved Bloom filter generation")
