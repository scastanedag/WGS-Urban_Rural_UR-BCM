#!/bin/bash

# Revisar que se hayan proporcionado las carpetas de entrada y salida
if [ "$#" -ne 2 ]; then
    echo "Uso: $0 <carpeta_input> <carpeta_output>"
    exit 1
fi

# Carpeta de entrada y salida
input_dir="$1"
output_dir="$2"

# Crear la carpeta de salida si no existe
mkdir -p "$output_dir"

# Iterar sobre todos los archivos forward (*.1.unaligned.fq.gz) en la carpeta de entrada
for forward_file in "$input_dir"/*.1.unaligned.fq.gz; do
    # Obtener el nombre base de la muestra (antes del ".1.unaligned.fq.gz")
    sample_name=$(basename "$forward_file" .1.unaligned.fq.gz)
    
    # Construir el nombre del archivo reverse correspondiente
    reverse_file="$input_dir/${sample_name}.2.unaligned.fq.gz"
    
    # Verificar si el archivo reverse existe
    if [ ! -f "$reverse_file" ]; then
        echo "Archivo reverse no encontrado para ${sample_name}, omitiendo..."
        continue
    fi

    # Construir el nombre del archivo de salida
    output_file="$output_dir/${sample_name}.merged.fq.gz"
    
    # Concatenar los archivos forward y reverse
    cat "$forward_file" "$reverse_file" > "$output_file"
    
    echo "Archivos ${forward_file} y ${reverse_file} concatenados en ${output_file}"
done


# para ejecutar:
# ./merge_samples_with_dirs.sh /volumen1/sergio/Metagenomics_Roj/J201006_NB552064_0188_AHNCN7BGXG_shep/clean_reads /volumen1/sergio/Metagenomics_Roj/J201006_NB552064_0188_AHNCN7BGXG_shep/clean_reads/merged
# ejecutar humann 
# for f in *.fq.gz; do humann -i $f -o humann_F; done
