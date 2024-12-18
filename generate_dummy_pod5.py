import pod5
import numpy as np
import random
import os
import uuid

def generate_random_sequence(length):
    """Generate a random DNA sequence"""
    bases = ['A', 'C', 'G', 'T']
    return ''.join(random.choice(bases) for _ in range(length))

def create_dummy_pod5_files(num_files=5, reads_per_file=100, read_length=1000):
    """
    Create dummy pod5 files with random sequences
    
    :param num_files: Number of pod5 files to generate
    :param reads_per_file: Number of reads in each file
    :param read_length: Length of each read
    """
    dummy_data_dir = 'dummy_data'
    os.makedirs(dummy_data_dir, exist_ok=True)
    
    for i in range(num_files):
        filename = os.path.join(dummy_data_dir, f'dummy_reads_{i}.pod5')
        
        with pod5.Writer(filename) as writer:
            for _ in range(reads_per_file):
                # Generate random read data
                sequence = generate_random_sequence(read_length)
                quality = np.random.randint(0, 40, size=read_length, dtype=np.uint8)
                signal = np.random.randint(0, 1000, size=read_length * 5, dtype=np.int16)
                
                # Create a read
                writer.add_read(
                    read_id=str(uuid.uuid4()),
                    signal=signal,
                    pore=1,
                    calibration=pod5.Calibration(
                        offset=0.0,
                        scale=1.0
                    ),
                    end_reason=1,
                    run_info=pod5.RunInfo(
                        acquisition_id=str(uuid.uuid4()),
                        acquisition_start_time=0.0,
                        adc_max=2048,
                        adc_min=0,
                        context_tags={},
                        experiment_name="dummy_experiment",
                        flow_cell_id="dummy_flow_cell",
                        flow_cell_product_code="dummy_product_code",
                        protocol_name="dummy_protocol",
                        protocol_run_id=str(uuid.uuid4()),
                        protocol_start_time=0.0,
                        sample_id="dummy_sample",
                        sample_rate=4000,
                        sequencing_kit="dummy_kit",
                        sequencer_position="dummy_position",
                        sequencer_position_type="dummy_type",
                        software="dummy_software",
                        system_name="dummy_system",
                        system_type="dummy_type",
                        tracking_id={},
                    ),
                    read_number=1,
                    start_sample=0,
                    median_before=0.0,
                    tracked_scaling=True,
                    predicted_scaling=True,
                    num_reads_since_mux_change=1,
                    time_since_mux_change=0.0,
                    sequence=sequence,
                    qstring=bytes([33 + q for q in quality])
                )
        
        print(f"Created {filename}")

if __name__ == "__main__":
    create_dummy_pod5_files()
