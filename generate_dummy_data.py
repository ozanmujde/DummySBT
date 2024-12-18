import os
import random
import json

def generate_random_sequence(length):
    """Generate a random DNA sequence"""
    bases = ['A', 'C', 'G', 'T']
    return ''.join(random.choice(bases) for _ in range(length))

def create_dummy_data(num_files=5, reads_per_file=100, read_length=1000):
    """
    Create dummy sequence files with random DNA sequences
    
    :param num_files: Number of files to generate
    :param reads_per_file: Number of reads in each file
    :param read_length: Length of each read
    """
    dummy_data_dir = 'dummy_data'
    os.makedirs(dummy_data_dir, exist_ok=True)
    
    for i in range(num_files):
        filename = os.path.join(dummy_data_dir, f'dummy_reads_{i}.json')
        
        reads = []
        for j in range(reads_per_file):
            read = {
                'id': f'read_{i}_{j}',
                'sequence': generate_random_sequence(read_length),
                'quality': [random.randint(0, 40) for _ in range(read_length)]
            }
            reads.append(read)
        
        with open(filename, 'w') as f:
            json.dump(reads, f)
        
        print(f"Created {filename}")

if __name__ == "__main__":
    create_dummy_data()
