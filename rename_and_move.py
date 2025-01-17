import os
import shutil
import argparse

def rename_and_move_files(input_dir, output_dir):
    # Verificar que la carpeta de entrada exista
    if not os.path.exists(input_dir):
        print(f"Error: El directorio de entrada '{input_dir}' no existe.")
        return

    # Crear la carpeta de salida si no existe
    os.makedirs(output_dir, exist_ok=True)

    # Procesar archivos en la carpeta de entrada
    for file_name in os.listdir(input_dir):
        if file_name.endswith(".1_paired_trim.fq.gz"):
            new_name = file_name.replace(".1_paired_trim.fq.gz", "_1.fq.gz")
            shutil.move(os.path.join(input_dir, file_name), os.path.join(output_dir, new_name))
        elif file_name.endswith(".2_paired_trim.fq.gz"):
            new_name = file_name.replace(".2_paired_trim.fq.gz", "_2.fq.gz")
            shutil.move(os.path.join(input_dir, file_name), os.path.join(output_dir, new_name))

    print(f"Archivos renombrados y movidos a '{output_dir}'.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Renombrar y mover archivos FASTQ.")
    parser.add_argument("-i", "--input", required=True, help="Directorio de entrada con los archivos.")
    parser.add_argument("-o", "--output", required=True, help="Directorio de salida donde se moverán los archivos.")
    args = parser.parse_args()

    rename_and_move_files(args.input, args.output)

#ejecucion
#python rename_and_move.py -i <ruta_a_la_carpeta_de_entrada> -o <ruta_a_la_carpeta_de_salida>
