import time
import psutil
import os
import shutil
import numpy as np
import json
import random
from memory_profiler import profile
from rawhash import RawHash
from sequence_bloom_tree import SequenceBloomTree
from generate_dummy_data import create_dummy_data

class SBTBenchmark:
    def __init__(self):
        self.process = psutil.Process()
        # Use fixed parameters for consistency
        bloom_size = 1048576  # 1MB Bloom filter
        hash_count = 3
        k_mer_size = 20
        
        self.rawhash = RawHash(k_mer_size=k_mer_size, bloom_size=bloom_size, hash_count=hash_count)
        self.sbt = SequenceBloomTree(k_mer_size=k_mer_size, bloom_size=bloom_size, hash_count=hash_count)
        
    def measure_memory(self):
        """Measure current memory usage in MB"""
        return self.process.memory_info().rss / 1024 / 1024
    
    def cleanup_dummy_data(self):
        """Remove existing dummy data"""
        if os.path.exists('dummy_data'):
            shutil.rmtree('dummy_data')
    
    def load_sequences_from_json(self, filepath):
        """Load sequences from a JSON file"""
        with open(filepath, 'r') as f:
            data = json.load(f)
            return [read['sequence'] for read in data]
    
    @profile
    def benchmark_sequence_processing(self, num_files=5, reads_per_file=100):
        """Benchmark sequence file processing"""
        print("\n=== Sequence Processing Benchmark ===")
        
        # Check if data already exists
        if not os.path.exists('dummy_data'):
            print("Generating new test data...")
            start_time = time.time()
            create_dummy_data(num_files, reads_per_file)
            data_gen_time = time.time() - start_time
            print(f"Data generation time: {data_gen_time:.2f} seconds")
        else:
            print("Using existing test data...")
        
        # Process files and measure performance
        start_mem = self.measure_memory()
        start_time = time.time()
        
        signatures = []
        sequences = []
        
        # Only process the number of files we need
        file_count = 0
        for filename in sorted(os.listdir('dummy_data')):
            if filename.endswith('.json') and file_count < num_files:
                filepath = os.path.join('dummy_data', filename)
                file_sequences = self.load_sequences_from_json(filepath)
                if len(file_sequences) >= reads_per_file:
                    # Only take the number of reads we need
                    file_sequences = file_sequences[:reads_per_file]
                    file_signatures = self.rawhash.process_file(filepath)[:reads_per_file]
                    signatures.extend(file_signatures)
                    sequences.extend(file_sequences)
                    file_count += 1
                    
        if file_count < num_files:
            print(f"Warning: Only found {file_count} valid files, needed {num_files}")
        
        processing_time = time.time() - start_time
        memory_used = self.measure_memory() - start_mem
        
        print(f"Processing time: {processing_time:.2f} seconds")
        print(f"Memory used: {memory_used:.2f} MB")
        print(f"Average time per read: {processing_time/(file_count*reads_per_file):.4f} seconds")
        
        return signatures, sequences
    
    @profile
    def benchmark_tree_construction(self, signatures, sequences):
        """Benchmark SBT construction"""
        print("\n=== Tree Construction Benchmark ===")
        
        start_mem = self.measure_memory()
        start_time = time.time()
        
        self.sbt.insert_signatures(signatures, sequences)
        
        construction_time = time.time() - start_time
        memory_used = self.measure_memory() - start_mem
        
        print(f"Construction time: {construction_time:.2f} seconds")
        print(f"Memory used: {memory_used:.2f} MB")
        print(f"Average time per signature: {construction_time/len(signatures):.4f} seconds")
    
    def mutate_signature(self, signature, num_mutations):
        """
        Create a mutated version of a Bloom filter signature
        
        :param signature: Original Bloom filter signature
        :param num_mutations: Number of bit flips to perform
        :return: Mutated signature
        """
        mutated = signature.copy()
        size = len(signature)
        # Randomly flip bits
        flip_positions = random.sample(range(size), k=num_mutations)
        for pos in flip_positions:
            mutated[pos] = not mutated[pos]
        return mutated
    
    @profile
    def benchmark_querying(self, signatures, sequences, num_queries=100):
        """Benchmark query performance"""
        print("\n=== Query Performance Benchmark ===")
        
        start_mem = self.measure_memory()
        start_time = time.time()
        
        # Initialize metrics
        metrics_sum = {
            'true_positives': 0,
            'false_positives': 0,
            'false_negatives': 0,
            'precision': 0,
            'recall': 0,
            'f1_score': 0
        }
        
        # Use existing SBT instance
        self.sbt.insert_signatures(signatures, sequences)
        
        # Generate and run queries
        for i in range(num_queries):
            # Randomly select a sequence to query
            idx = random.randint(0, len(signatures) - 1)
            query_sig = signatures[idx]
            
            # Add some mutations if it's not an exact match query
            if i >= num_queries // 2:  # 50% mutated queries
                num_mutations = random.randint(1, 5)
                query_sig = self.mutate_signature(query_sig, num_mutations)
            
            # Run query with expected matches
            expected = [idx]  # The original sequence should match
            result = self.sbt.query(query_sig, expected_matches=expected)
            
            # Update metrics
            for key in metrics_sum:
                metrics_sum[key] += result['metrics'][key]
        
        # Calculate average metrics
        for key in metrics_sum:
            metrics_sum[key] /= num_queries
        
        query_time = time.time() - start_time
        memory_used = self.measure_memory() - start_mem
        
        print(f"\nQuery Performance Metrics:")
        print(f"Average precision: {metrics_sum['precision']:.3f}")
        print(f"Average recall: {metrics_sum['recall']:.3f}")
        print(f"Average F1 score: {metrics_sum['f1_score']:.3f}")
        print(f"True positives: {metrics_sum['true_positives']:.1f}")
        print(f"False positives: {metrics_sum['false_positives']:.1f}")
        print(f"False negatives: {metrics_sum['false_negatives']:.1f}")
        print(f"\nPerformance Metrics:")
        print(f"Total query time: {query_time:.2f} seconds")
        print(f"Average time per query: {query_time/num_queries:.4f} seconds")
        print(f"Memory used: {memory_used:.2f} MB")
    
    def run_full_benchmark(self, num_files=5, reads_per_file=100, num_queries=100):
        """Run complete benchmark suite"""
        print("\n=== Starting Full Benchmark ===")
        print(f"Configuration:")
        print(f"- Number of files: {num_files}")
        print(f"- Reads per file: {reads_per_file}")
        print(f"- Number of queries: {num_queries}")
        
        # Store results in a file for comparison
        result_file = "benchmark_results.txt"
        with open(result_file, "a") as f:
            f.write(f"\n=== Benchmark Run at {time.strftime('%Y-%m-%d %H:%M:%S')} ===\n")
            f.write(f"Files: {num_files}, Reads/File: {reads_per_file}, Queries: {num_queries}\n")
        
        signatures, sequences = self.benchmark_sequence_processing(num_files, reads_per_file)
        self.benchmark_tree_construction(signatures, sequences)
        self.benchmark_querying(signatures, sequences, num_queries)
        
        print("\n=== Benchmark Complete ===")
        print(f"Results saved to {result_file}")

if __name__ == "__main__":
    # Add memory-profiler and psutil to requirements
    with open("requirements.txt", "a") as f:
        f.write("\npsutil==5.9.8\nmemory-profiler==0.61.0\n")
    
    # Run benchmarks with different configurations
    benchmark = SBTBenchmark()
    
    # # Small dataset
    # print("\n=== Small Dataset Benchmark ===")
    # benchmark.run_full_benchmark(num_files=2, reads_per_file=50, num_queries=50)
    
    # # Medium dataset
    # print("\n=== Medium Dataset Benchmark ===")
    # benchmark.run_full_benchmark(num_files=5, reads_per_file=100, num_queries=100)
    
    # Large dataset
    print("\n=== Large Dataset Benchmark ===")
    benchmark.run_full_benchmark(num_files=10, reads_per_file=200, num_queries=200)
