##!/bin/sh

OUT_DIR=$1
REF_GENOME=$2
OUT_FILE=$3
FASTP_OUT=$(ls $OUT_DIR*json)
BAM_FILE=$(ls $OUT_DIR*sorted.bam)
SNPEFF_FILE=$(ls $OUT_DIR*ann.vcf)

echo "" >> $OUT_FILE
echo -e "Respostas do Desafio Técnico 5\n" > $OUT_FILE

echo "1 - Quantas sequências de DNA de paciente sequenciados temos nos arquivos de fastqs R1 e R2 respectivamente ?" >> $OUT_FILE
R1_READS=$(grep -A 1 "read1_before_filtering" $FASTP_OUT | tail -n 1 | sed -e 's/.*: //g' -e 's/,//g')
R2_READS=$(grep -A 1 "read2_before_filtering" $FASTP_OUT | tail -n 1 | sed -e 's/.*: //g' -e 's/,//g')
echo -e "R: Existem $R1_READS sequencias no arquivo R1 e $R2_READS sequencias no arquivo R2. Totalizando $(($R1_READS + $R2_READS)) sequencias.\n" >> $OUT_FILE

echo "2 - Sobre o genoma humano hg19, quantos contigs tem o nosso genoma hg19 (hg19.fasta) aqui disponibilizado para este pipeline ?" >> $OUT_FILE
echo -e "R: Existem $(grep -c ">" $REF_GENOME) contigs no arquivo $REF_GENOME.\n" >> $OUT_FILE

echo "3 - Quantos alinhamentos há na região chr17:41197694-41197819 ?" >> $OUT_FILE
echo -e "chr17\t41197694\t41197819" > $OUT_DIR/ch17.bed
echo -e "R: Existem $(samtools view -c -L $OUT_DIR/ch17.bed $BAM_FILE) sequencias mapeadas na região chr17:41197694-41197819\n" >> $OUT_FILE

echo "4 - Quantos alinhamentos não conseguiram ser mapeados (unmapped alignments ?)" >> $OUT_FILE
UNALGN_TOTAL=$(samtools view -f 4 -c $BAM_FILE)
echo -e "R: $UNALGN_TOTAL sequencias não foram mapeadas considerando todo o genoma de referência.\n" >> $OUT_FILE

echo -e "5 - Realize o alinhamento das sequências FASTQ contra o genoma de referência hg19 aqui disponibilizado, e realize a chamada de variantes utilizando a região alvo BRCA.list (interesse apenas na região dos genes BRCA1 e BRCA2) Realize a anotação de dados funcionais usando o SNPEFF.Com este arquivo em mãos, responda as seguintes perguntas:\n" >> $OUT_FILE
#criando uma tabela com o sumário das informações
echo -e "chromossome\tposition\tref\talt\timpact" > $OUT_DIR/impact_table.tsv
#Apesar de uma mesma variante ter poder ter mais de um impacto anotado, todos os impactos estão na mesma linha. Analisando pelo menos estes dados, não existe uma mesma variante com 2 tipos de impacto diferente, por exemplo, efeito HIGH e MODERATE na mesma variante, então optei por uma busca simples com o grep.
grep "HIGH" $SNPEFF_FILE | cut -f 1,2,4,5 | awk -v OFS='\t' '{print $0,"HIGH"}' >> $OUT_DIR/impact_table.tsv
grep "MODERATE" $SNPEFF_FILE | cut -f 1,2,4,5 | awk -v OFS='\t' '{print $0,"MODERATE"}' >> $OUT_DIR/impact_table.tsv
grep "LOW" $SNPEFF_FILE | cut -f 1,2,4,5 | awk -v OFS='\t' '{print $0,"LOW"}' >> $OUT_DIR/impact_table.tsv
echo "5.1 - Quantas variantes existem com impacto funcional (Annotation_Impact) do tipo HIGH, MODERATE, LOW ?" >> $OUT_FILE
echo -e "R: Existem $(grep -c 'HIGH' $OUT_DIR/impact_table.tsv) variante do tipo HIGH, $(grep -c "MODERATE" $OUT_DIR/impact_table.tsv) variantes do tipo MODERATE, e $(grep -c "LOW" $OUT_DIR/impact_table.tsv) variantes do tipo LOW.\n" >> $OUT_FILE
echo "5.2 - Existe alguma variante em HIGH ? Qual é cromossomo, posição e a alteração ?" >> $OUT_FILE
echo "R: A variante do tipo HIGH se encontra no cromossomo $(grep "HIGH" $OUT_DIR/impact_table.tsv | cut -f 1), na posição $(grep "HIGH" $OUT_DIR/impact_table.tsv | cut -f 2), onde ocorreu uma alteração de $(grep "HIGH" $OUT_DIR/impact_table.tsv | cut -f 3) para $(grep "HIGH" $OUT_DIR/impact_table.tsv | cut -f 4)." >> $OUT_FILE
