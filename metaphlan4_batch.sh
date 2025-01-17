#!/bin/bash

#se instala metaphlan4 siguiendo https://github.com/biobakery/MetaPhlAn/wiki/MetaPhlAn-4.1#installation
#la base se instala en /bases

#para usar 
#metaphlan /path/to/sample_R1.fastq.gz,/path/to/sample_R2.fastq.gz \
#    --input_type fastq \
#    --bowtie2db /path/to/metaphlan_db \
#    --output_file /path/to/output/profiled_metagenome.txt \
#    --bowtie2out /path/to/output/bowtie2_output.bz2 \
#    --nproc 8


# Ruta al directorio con las muestras de reads no alineados. Ajustar de acuerdo con la ruta
INPUT_DIR="/volumen1/sergio/Metagenomics_Roj/F200921_NS500355_0991_AHGVM2BGXG_shep/clean_reads"
# Ruta a la base de datos de MetaPhlAn4. Ajustar de acuerdo con la ruta
METAPHLAN_DB="/volumen1/sergio/Metagenomics_Roj/bases"
# Ruta al directorio de salida. Ajustar de acuerdo con la ruta
OUTPUT_DIR="/volumen1/sergio/Metagenomics_Roj/F200921_NS500355_0991_AHGVM2BGXG_shep/clean_reads/Metaphlan"
# Número de hilos a usar
THREADS=8

# Crear el directorio de salida si no existe
mkdir -p "$OUTPUT_DIR"

# Procesar todas las muestras
for R1_FILE in "$INPUT_DIR"/*.1.unaligned.fq.gz; do
    # Derivar el nombre base de la muestra
    SAMPLE_NAME=$(basename "$R1_FILE" .1.unaligned.fq.gz)
    
    # Derivar la ruta al archivo R2 correspondiente
    R2_FILE="${INPUT_DIR}/${SAMPLE_NAME}.2.unaligned.fq.gz"
    
    # Verificar que ambos archivos existan
    if [[ -f "$R1_FILE" && -f "$R2_FILE" ]]; then
        echo "Procesando muestra: $SAMPLE_NAME"
        
        # Definir los archivos de salida
        OUTPUT_PROFILE="${OUTPUT_DIR}/${SAMPLE_NAME}_profiled_metagenome.txt"
        OUTPUT_BOWTIE2="${OUTPUT_DIR}/${SAMPLE_NAME}_bowtie2_output.bz2"
        
        # Ejecutar MetaPhlAn4
        metaphlan "$R1_FILE","$R2_FILE" \
            --input_type fastq \
            --bowtie2db "$METAPHLAN_DB" \
            --output_file "$OUTPUT_PROFILE" \
            --bowtie2out "$OUTPUT_BOWTIE2" \
            --nproc "$THREADS"
        
        echo "Muestra $SAMPLE_NAME procesada."
    else
        echo "Archivos R1 o R2 no encontrados para la muestra: $SAMPLE_NAME"
    fi
done

echo "Procesamiento completado."
