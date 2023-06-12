
Example data for CSSRLIB

# QZSS CLAS and SLAS

The correction data for CLAS and SLAS is downloaded from QZS Web site (https://qzss.go.jp/en/).
Please check the following license term for these data set.

[correction data for CLAS]
https://sys.qzss.go.jp/dod/en/archives/agree.html?SERVICE_ID=CLAS

[correction data for SLAS]
https://sys.qzss.go.jp/dod/en/archives/agree.html?SERVICE_ID=SLAS

# Galileo HAS

The files for the Galileo-HAS test script Galileo-HAS-SIS-ICD-1.0_Annex_D_HAS_Message_Decoding_Example.txt and Galileo-HAS-SIS-ICD_1.0_Annex_B_Reed_Solomon_Generator_Matrix.txt are embedded as appendix in the Galileo HAS-ICD.

The file with the HAS message decoding example in Annex D contains non-UTF-8 characters. These must be removed with

```
iconv -f utf-8 -t utf-8 -c Galileo-HAS-SIS-ICD-1.0_Annex_D_HAS_Message_Decoding_Example.txt > Galileo-HAS-SIS-ICD-1.0_Annex_D_HAS_Message_Decoding_Example.fixed.txt
```