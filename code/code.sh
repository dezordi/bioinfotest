##!/bin/sh

: ' Não vou inserir a etapa de indexação do genoma de referencia, pois nos dados de download os arquivos de indexação ja estão presentes.
Em outros casos, seria necessário rodar: bwa index <REF_GENOME> '

#criando argumentos
REF_GENOME=$1        #genoma de referencia
BED_LIST=$2          #arquivo bed com regioes de interesse
FQ_1=$3              #arquivo fastq R1
FQ_2=$4              #arquivo fastq R2
THREADS=$5           #número de threads
OUT_DIR=$6           #diretório de output

#criando prefixo para os arquivos de saída
PREFIX=$(echo $FQ_1 | sed -e 's/.*\///g'| awk -F"_L001" '{print $1}') #sed para remover o caminho absoluto ou relativo, awk para pegar o codigo da amostra.

: ' Apesar de não estar nas análises solicitadas, vou inserir uma etapa de tratamento de dados.
    Vou utilizar a ferramenta fastp por diversos fatores:
        - Engloba a etapa de analise de qualidade e trimmagem dos dados;
        - É mais rapida que as outras soluções;
        - Detecta automaticamente adaptadores e os remove.
        - Utilizei o parametro de qualidade média minima na janela de análise para 30, que em PHRED score representa 99,9% de acurácia na base sequênciada'

fastp -i $FQ_1 -I $FQ_2 -w $THREADS -o $OUT_DIR$PREFIX.R1.fq.gz -O $OUT_DIR$PREFIX.R2.fq.gz -h $OUT_DIR$PREFIX.fastp.html -j $OUT_DIR$PREFIX.fastp.json --cut_front --qualified_quality_phred 30

bwa mem $REF_GENOME $OUT_DIR$PREFIX.R1.fq.gz $OUT_DIR$PREFIX.R2.fq.gz > $OUT_DIR$PREFIX.bam

samtools sort $OUT_DIR$PREFIX.bam > $OUT_DIR$PREFIX.sorted.bam
samtools index $OUT_DIR$PREFIX.sorted.bam

freebayes -f $REF_GENOME $OUT_DIR$PREFIX.sorted.bam --targets $BED_LIST > $OUT_DIR$PREFIX.vcf

#Foi passado o argumento -Xmx4g, pois no default poderia dar limite de memória no java dependente do ambiente de execução.
snpEff ann -Xmx4g -noStats hg19 $OUT_DIR$PREFIX.vcf > $OUT_DIR$PREFIX.ann.vcf
