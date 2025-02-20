import time

def write(initial_data, results, output_filename):
    start_time = time.time()

    with open(output_filename, 'w') as output_file:
        # Header message
        header_message = """Photonics UdeMedellín \n
Version 1 - 2024 \n

Calculations for photonic crystals in Python \n
Date: {}
\n""".format(time.strftime("%Y-%m-%d %H:%M:%S"))
        output_file.write(header_message)

        output_file.write("Initial data of the problem:\n")
        output_file.write(f"{initial_data}\n\n")

        output_file.write("Results obtained:\n")
        for result in results:
            output_file.write(f"{result}\n")
        output_file.write("\n")

    end_time = time.time()
    total_time = end_time - start_time

    with open(output_filename, 'a') as output_file:  # 'a' to open in 'append' mode
        output_file.write(f"Total calculation time: {total_time:.4f} seconds\n")
