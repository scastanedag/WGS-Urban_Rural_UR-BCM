
#!/bin/bash

# Carpeta donde se encuentran las muestras
input_folder="/ruta/a/la/carpeta/de/las/muestras"

# Carpeta donde se guardarán los resultados
output_folder="/ruta/a/la/carpeta/de/resultados"

# Entra al directorio de entrada
cd "$input_folder" || exit

# Itera sobre los archivos .1.fastq.gz en el directorio
for file in *R1.fastq.gz; do
    # Extrae el nombre de la muestra (sin la extensión)
    sample_name=$(basename "$file" _R1.fastq.gz)

    # Define los nombres de los archivos de entrada para FastQC
    input_file1="$input_folder/$sample_name"_R1.fastq.gz
    input_file2="$input_folder/$sample_name"_R2.fastq.gz

    # Define el nombre del directorio de salida para FastQC
    output_dir="$output_folder/$sample_name"_fastqc

    # Crea el directorio de salida si no existe
    mkdir -p "$output_dir"

    # Ejecuta FastQC para las muestras de lectura 1 y lectura 2
    fastqc -o "$output_dir" "$input_file1" "$input_file2"
done
