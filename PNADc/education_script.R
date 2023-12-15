library('PNADcIBGE')
library('survey')

pnad_2022 <- get_pnadc(year=2022, interview=5)

filtered_design <- subset(pnad_2022, V2009 >= 16)

total_education_1 <- svytotal(x=~interaction(V3003A, V1023, UF),
                            design = filtered_design,
                            na.rm = TRUE)

write.csv(total_education_1,
          "total_education_1_V1023.csv",
          row.names=TRUE)

total_education_2 <- svytotal(x=~interaction(V3009A, V1023, UF),
                              design = filtered_design,
                              na.rm = TRUE)

write.csv(total_education_2,
          "total_education_2_V1023.csv",
          row.names=TRUE)

filtered_design <- subset(filtered_design, is.na(V3003A))
filtered_design <- subset(filtered_design, is.na(V3009A))
filtered_design <- subset(filtered_design, !is.na(V3001))

total_education_3 <- svytotal(x=~interaction(V3001, V1023, UF),
                              design = filtered_design,
                              na.rm = TRUE)

write.csv(total_education_3,
          "total_education_3_V1023.csv",
          row.names=TRUE)

print(total_education_1)
print(total_education_2)
print(total_education_3)