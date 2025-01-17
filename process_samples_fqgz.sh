#!/bin/bash
# Script para alinear lecturas contra el genoma humano y extraer las lecturas no alineadas para múltiples muestras en formato .fastq.gz.
# Ahora permite ingresar las rutas del índice y el genoma como parámetros.

set -e

# Comprobar que se ingresaron los parámetros necesarios
if [ $# -lt 3 ]; then
    echo "Uso: $0 <directorio_entrada> <directorio_salida> <ruta_genoma_humano> <ruta_indice>"
    echo "Ejemplo: $0 /ruta/a/lecturas /ruta/salida /ruta/a/genoma_humano.fasta /ruta/a/indice_genoma"
    exit 1
fi

# Directorios de entrada y salida
input_dir="$1"  # El primer argumento es el directorio de entrada que contiene las muestras .fastq.gz
output_dir="$2"  # El segundo argumento es el directorio de salida donde se guardarán los resultados

# Ruta al genoma humano y el índice de Bowtie2 proporcionados como argumentos
human_genome="$3"  # Genoma humano (archivo FASTA)
human_genome_index="$4"  # Índice de Bowtie2 (sin extensión, solo el prefijo)

# Verificar si el índice existe, si no, lo crea
if [ ! -f "${human_genome_index}.1.bt2" ]; then
    echo "El índice de Bowtie2 no existe. Creando el índice para el genoma humano..."
    bowtie2-build "$human_genome" "$human_genome_index"
fi

# Crear el directorio de salida si no existe
mkdir -p "$output_dir"

# Procesar todas las muestras .fastq.gz en el directorio de entrada
for read_file in "$input_dir"/*.fastq.gz; do
    # Comprobar que el archivo es un archivo .fastq.gz
    if [[ "$read_file" == *.fastq.gz ]]; then
        echo "Procesando archivo: $read_file"
        
        # Obtener el nombre base del archivo (sin la extensión .fastq.gz)
        sample_name=$(basename "$read_file" .fastq.gz)
        
        # Alinear las lecturas con el genoma humano
        echo "Alineando $sample_name con el genoma humano..."
        bowtie2 -x "$human_genome_index" -U "$read_file" -S "$output_dir/$sample_name.sam"
        
        # Filtrar las lecturas no alineadas
        echo "Filtrando lecturas no alineadas para $sample_name..."
        samtools view -b -f 4 "$output_dir/$sample_name.sam" > "$output_dir/$sample_name.unaligned.bam"
        
        # Convertir el archivo BAM a FASTQ
        echo "Convirtiendo BAM a FASTQ para $sample_name..."
        samtools fastq "$output_dir/$sample_name.unaligned.bam" | gzip > "$output_dir/$sample_name.unaligned.fastq.gz"
        
        # Limpiar archivos intermedios
        rm "$output_dir/$sample_name.sam"
        rm "$output_dir/$sample_name.unaligned.bam"
        
        echo "Muestra $sample_name procesada con éxito."
    else
        echo "El archivo $read_file no es un archivo .fastq.gz, se omite."
    fi
done

echo "Todos los archivos han sido procesados."


#./process_samples_fqgz.sh /volumen1/sergio/Metagenomics_Roj/F200921_NS500355_0991_AHGVM2BGXG_shep /volumen1/sergio/Metagenomics_Roj/F200921_NS500355_0991_AHGVM2BGXG_shep/clean_reads /volumen1/sergio/Metagenomics_Roj/Mixed.fasta /volumen1/sergio/Metagenomics_Roj/ref_index
