import os
from rawhash import RawHash
from sequence_bloom_tree import SequenceBloomTree
from generate_dummy_pod5 import create_dummy_pod5_files

def main():
    # Generate dummy pod5 files if they don't exist
    dummy_data_dir = 'dummy_data'
    if not os.path.exists(dummy_data_dir) or len(os.listdir(dummy_data_dir)) == 0:
        create_dummy_pod5_files()
    
    # Initialize RawHash and SequenceBloomTree
    rawhash = RawHash()
    sbt = SequenceBloomTree()
    
    # Process all dummy pod5 files
    for filename in os.listdir(dummy_data_dir):
        if filename.endswith('.pod5'):
            filepath = os.path.join(dummy_data_dir, filename)
            
            # Generate signatures for the file
            signatures = rawhash.process_pod5_file(filepath)
            
            # Insert signatures into Sequence Bloom Tree
            sbt.insert_signatures(signatures)
    
    # Demonstrate querying
    print("Sequence Bloom Tree Processing Complete!")
    print("Total signatures processed:", len(sbt.tree[0]))
    
    # Query with a random signature from the first file
    first_file = os.path.join(dummy_data_dir, os.listdir(dummy_data_dir)[0])
    query_signatures = rawhash.process_pod5_file(first_file)
    
    if query_signatures:
        query_result = sbt.query(query_signatures[0])
        print("Query Results:", query_result)

if __name__ == "__main__":
    main()
