library("PNADcIBGE")
library("survey")

df <- get_pnadc(year=2022, interview=5, deflator=2022, design=FALSE)

df <- subset(df, V2009 >= 18) # Seleciona apenas os adultos

df$soma_renda <- ave(df$VD4022, df$ID_DOMICILIO, FUN=function(x) sum(x, na.rm=TRUE))

#df$cat_renda <- cut(df$soma_renda, breaks = c(-Inf, 105, 210, 303, 909, 1212, 2424, 3636, 6060, 12120, Inf)) # Categorias em termos de S.M.
df$cat_renda <- cut(df$soma_renda, breaks = c(-Inf, 2000, 3000, 5000, 10000, Inf)) # Categorias em reais

df <- df[!duplicated(df$ID_DOMICILIO), ] # Conta o número de domicílios ao invés de pessoas

svy <- pnadc_design(df)
test <- svyby(formula=~cat_renda, by=~interaction(UF, V1023), FUN=svytotal, design=svy, na.rm=TRUE)
write.csv(test, "income_Capital_npeople_minwage.csv", row.names=TRUE)