## Making Complete Randomised Complete Block Design for experiment layout
# 10/02/2020 

library(agricolae)

# This works
#It is treating it as 9 rows of 9, I think just randomising within a row. So could never see two
# of a given species in the same row
trt <- c("VERO", "ARCA", "PLDE", "POLE", "PEAI", "HYGL", "TROR", "TRCY", "LARO")
outdesign <- design.rcbd(trt, 9, serie=2, 1234, "default") #seed is fourth argument here
# Not sure what serie does or how randomisation types differ (what "default" is)
print(outdesign$sketch)


#Problem with below is it is treating it as 3 rows of 27
trt <- c("VERO C", "ARCA C", "PLDE C", "POLE C", "PEAI C", "HYGL C", "TROR C", "TRCY C", "LARO C",
         "VERO E", "ARCA E", "PLDE E", "POLE E", "PEAI E", "HYGL E", "TROR E", "TRCY E", "LARO E",
         "VERO T", "ARCA T", "PLDE T", "POLE T", "PEAI T", "HYGL T", "TROR T", "TRCY T", "LARO T")
outdesign <- design.rcbd(trt, 3, serie=2, 1234, "default") #seed is fourth argument here
book <- outdesign$book 
write.table(book, "rcbd.txt", row.names=TRUE, sep="\t")
file.show("rcbd.txt")
fieldbook <- zigzag(outdesign)
print(outdesign$sketch)
print(matrix(fieldbook[,1], byrow=TRUE, ncol=9))
print(matrix(outdesign$sketch), byrow=TRUE, ncol=9)

outdesign <- designRandomize

#(trt, 3, serie=2, 1234, "default") #seed is fourth argument here

design.crd()

r <- c(3, 3, 3, 3, 3, 3, 3, 3, 3)
outdesign <- design.crd(trt, r, serie=2, 1234, "default") #seed is fourth argument here
book <- outdesign$book
print(outdesign$parameters)
write.table(book, "crd.txt")




#### With Aubrie in Carnamah, 14/02/20
##Setting seed guide - sites 1-8, A = 1, B = 2, C = 3. e.g. site 1 control (C) = 13.
## e.g. site 4 'B' treatment = 42

Subplots <- c("VERO_C1", "ARCA_C1", "PLDE_C1", "POLE_C1", "PEAI_C1", "HYGL_C1", "TROR_C1", "TRCY_C1", "LARO_C1",
              "VERO_C2", "ARCA_C2", "PLDE_C2", "POLE_C2", "PEAI_C2", "HYGL_C2", "TROR_C2", "TRCY_C2", "LARO_C2",
              "VERO_C3", "ARCA_C3", "PLDE_C3", "POLE_C3", "PEAI_C3", "HYGL_C3", "TROR_C3", "TRCY_C3", "LARO_C3",
              "VERO_E1", "ARCA_E1", "PLDE_E1", "POLE_E1", "PEAI_E1", "HYGL_E1", "TROR_E1", "TRCY_E1", "LARO_E1",
              "VERO_E2", "ARCA_E2", "PLDE_E2", "POLE_E2", "PEAI_E2", "HYGL_E2", "TROR_E2", "TRCY_E2", "LARO_E2",
              "VERO_E3", "ARCA_E3", "PLDE_E3", "POLE_E3", "PEAI_E3", "HYGL_E3", "TROR_E3", "TRCY_E3", "LARO_E3",
              "VERO_T1", "ARCA_T1", "PLDE_T1", "POLE_T1", "PEAI_T1", "HYGL_T1", "TROR_T1", "TRCY_T1", "LARO_T1",
              "VERO_T2", "ARCA_T2", "PLDE_T2", "POLE_T2", "PEAI_T2", "HYGL_T2", "TROR_T2", "TRCY_T2", "LARO_T2")
set.seed(11)
RandomisedVector <- sample(Subplots, 72, replace = FALSE)
MatrixVector <- matrix(RandomisedVector, nrow = 8, ncol = 9)
#View(MatrixVector)
write.csv(MatrixVector, file = "Output/Site_layouts_2020/Site_1_A.csv")
set.seed(12)
RandomisedVector <- sample(Subplots, 72, replace = FALSE)
MatrixVector <- matrix(RandomisedVector, nrow = 8, ncol = 9)
#View(MatrixVector)
write.csv(MatrixVector, file = "Output/Site_layouts_2020/Site_1_B.csv")
set.seed(13)
RandomisedVector <- sample(Subplots, 72, replace = FALSE)
MatrixVector <- matrix(RandomisedVector, nrow = 8, ncol = 9)
#View(MatrixVector)
write.csv(MatrixVector, file = "Output/Site_layouts_2020/Site_1_C.csv")
set.seed(21)
RandomisedVector <- sample(Subplots, 72, replace = FALSE)
MatrixVector <- matrix(RandomisedVector, nrow = 8, ncol = 9)
#View(MatrixVector)
write.csv(MatrixVector, file = "Output/Site_layouts_2020/Site_2_A.csv")

### Accidentally made a mistake so the seed for Site 2 B is 32
set.seed(32)
RandomisedVector <- sample(Subplots, 72, replace = FALSE)
MatrixVector <- matrix(RandomisedVector, nrow = 8, ncol = 9)
#View(MatrixVector)
write.csv(MatrixVector, file = "Output/Site_layouts_2020/Site_2_B.csv")
set.seed(23)
RandomisedVector <- sample(Subplots, 72, replace = FALSE)
MatrixVector <- matrix(RandomisedVector, nrow = 8, ncol = 9)
#View(MatrixVector)
write.csv(MatrixVector, file = "Output/Site_layouts_2020/Site_2_C.csv")
set.seed(31)
RandomisedVector <- sample(Subplots, 72, replace = FALSE)
MatrixVector <- matrix(RandomisedVector, nrow = 8, ncol = 9)
#View(MatrixVector)
write.csv(MatrixVector, file = "Output/Site_layouts_2020/Site_3_A.csv")

### Made a mistake so the seed for Site 3 B is 132
set.seed(132)
RandomisedVector <- sample(Subplots, 72, replace = FALSE)
MatrixVector <- matrix(RandomisedVector, nrow = 8, ncol = 9)
#View(MatrixVector)
write.csv(MatrixVector, file = "Output/Site_layouts_2020/Site_3_B.csv")
set.seed(33)
RandomisedVector <- sample(Subplots, 72, replace = FALSE)
MatrixVector <- matrix(RandomisedVector, nrow = 8, ncol = 9)
#View(MatrixVector)
write.csv(MatrixVector, file = "Output/Site_layouts_2020/Site_3_C.csv")
set.seed(41)
RandomisedVector <- sample(Subplots, 72, replace = FALSE)
MatrixVector <- matrix(RandomisedVector, nrow = 8, ncol = 9)
#View(MatrixVector)
write.csv(MatrixVector, file = "Output/Site_layouts_2020/Site_4_A.csv")
set.seed(42)
RandomisedVector <- sample(Subplots, 72, replace = FALSE)
MatrixVector <- matrix(RandomisedVector, nrow = 8, ncol = 9)
#View(MatrixVector)
write.csv(MatrixVector, file = "Output/Site_layouts_2020/Site_4_B.csv")
set.seed(43)
RandomisedVector <- sample(Subplots, 72, replace = FALSE)
MatrixVector <- matrix(RandomisedVector, nrow = 8, ncol = 9)
#View(MatrixVector)
write.csv(MatrixVector, file = "Output/Site_layouts_2020/Site_4_C.csv")
set.seed(51)
RandomisedVector <- sample(Subplots, 72, replace = FALSE)
MatrixVector <- matrix(RandomisedVector, nrow = 8, ncol = 9)
#View(MatrixVector)
write.csv(MatrixVector, file = "Output/Site_layouts_2020/Site_5_A.csv")
set.seed(52)
RandomisedVector <- sample(Subplots, 72, replace = FALSE)
MatrixVector <- matrix(RandomisedVector, nrow = 8, ncol = 9)
#View(MatrixVector)
set.seed(53)
RandomisedVector <- sample(Subplots, 72, replace = FALSE)
MatrixVector <- matrix(RandomisedVector, nrow = 8, ncol = 9)
#View(MatrixVector)
write.csv(MatrixVector, file = "Output/Site_layouts_2020/Site_5_C.csv")
set.seed(61)
RandomisedVector <- sample(Subplots, 72, replace = FALSE)
MatrixVector <- matrix(RandomisedVector, nrow = 8, ncol = 9)
#View(MatrixVector)
write.csv(MatrixVector, file = "Output/Site_layouts_2020/Site_6_A.csv")

set.seed(62)
RandomisedVector <- sample(Subplots, 72, replace = FALSE)
MatrixVector <- matrix(RandomisedVector, nrow = 8, ncol = 9)
#View(MatrixVector)
write.csv(MatrixVector, file = "Output/Site_layouts_2020/Site_6_B.csv")

set.seed(63)
RandomisedVector <- sample(Subplots, 72, replace = FALSE)
MatrixVector <- matrix(RandomisedVector, nrow = 8, ncol = 9)
#View(MatrixVector)
write.csv(MatrixVector, file = "Output/Site_layouts_2020/Site_6_C.csv")
set.seed(71)
RandomisedVector <- sample(Subplots, 72, replace = FALSE)
MatrixVector <- matrix(RandomisedVector, nrow = 8, ncol = 9)
#View(MatrixVector)
write.csv(MatrixVector, file = "Output/Site_layouts_2020/Site_7_A.csv")
set.seed(72)
RandomisedVector <- sample(Subplots, 72, replace = FALSE)
MatrixVector <- matrix(RandomisedVector, nrow = 8, ncol = 9)
#View(MatrixVector)
write.csv(MatrixVector, file = "Output/Site_layouts_2020/Site_7_B.csv")

Subplots <- c("VERO_C1", "ARCA_C1", "PLDE_C1", "POLE_C1", "PEAI_C1", "HYGL_C1", "TROR_C1", "TRCY_C1", "LARO_C1",
              "VERO_C2", "ARCA_C2", "PLDE_C2", "POLE_C2", "PEAI_C2", "HYGL_C2", "TROR_C2", "TRCY_C2", "LARO_C2",
              "VERO_C3", "ARCA_C3", "PLDE_C3", "POLE_C3", "PEAI_C3", "HYGL_C3", "TROR_C3", "TRCY_C3", "LARO_C3",
              "VERO_E1", "ARCA_E1", "PLDE_E1", "POLE_E1", "PEAI_E1", "HYGL_E1", "TROR_E1", "TRCY_E1", "LARO_E1",
              "VERO_E2", "ARCA_E2", "PLDE_E2", "POLE_E2", "PEAI_E2", "HYGL_E2", "TROR_E2", "TRCY_E2", "LARO_E2",
              "VERO_E3", "ARCA_E3", "PLDE_E3", "POLE_E3", "PEAI_E3", "HYGL_E3", "TROR_E3", "TRCY_E3", "LARO_E3",
              "VERO_T1", "ARCA_T1", "PLDE_T1", "POLE_T1", "PEAI_T1", "HYGL_T1", "TROR_T1", "TRCY_T1", "LARO_T1",
              "VERO_T2", "ARCA_T2", "PLDE_T2", "POLE_T2", "PEAI_T2", "HYGL_T2", "TROR_T2", "TRCY_T2", "LARO_T2")
set.seed(73)
RandomisedVector <- sample(Subplots, 72, replace = FALSE)
MatrixVector <- matrix(RandomisedVector, nrow = 8, ncol = 9)
#View(MatrixVector)
write.csv(MatrixVector, file = "Output/Site_layouts_2020/Site_7_C.csv")

Subplots <- c("VERO_C1", "ARCA_C1", "PLDE_C1", "POLE_C1", "PEAI_C1", "HYGL_C1", "TROR_C1", "TRCY_C1", "LARO_C1",
              "VERO_C2", "ARCA_C2", "PLDE_C2", "POLE_C2", "PEAI_C2", "HYGL_C2", "TROR_C2", "TRCY_C2", "LARO_C2",
              "VERO_C3", "ARCA_C3", "PLDE_C3", "POLE_C3", "PEAI_C3", "HYGL_C3", "TROR_C3", "TRCY_C3", "LARO_C3",
              "VERO_E1", "ARCA_E1", "PLDE_E1", "POLE_E1", "PEAI_E1", "HYGL_E1", "TROR_E1", "TRCY_E1", "LARO_E1",
              "VERO_E2", "ARCA_E2", "PLDE_E2", "POLE_E2", "PEAI_E2", "HYGL_E2", "TROR_E2", "TRCY_E2", "LARO_E2",
              "VERO_E3", "ARCA_E3", "PLDE_E3", "POLE_E3", "PEAI_E3", "HYGL_E3", "TROR_E3", "TRCY_E3", "LARO_E3",
              "VERO_T1", "ARCA_T1", "PLDE_T1", "POLE_T1", "PEAI_T1", "HYGL_T1", "TROR_T1", "TRCY_T1", "LARO_T1",
              "VERO_T2", "ARCA_T2", "PLDE_T2", "POLE_T2", "PEAI_T2", "HYGL_T2", "TROR_T2", "TRCY_T2", "LARO_T2")
set.seed(81)
RandomisedVector <- sample(Subplots, 72, replace = FALSE)
MatrixVector <- matrix(RandomisedVector, nrow = 8, ncol = 9)
#View(MatrixVector)
write.csv(MatrixVector, file = "Output/Site_layouts_2020/Site_8_A.csv")

Subplots <- c("VERO_C1", "ARCA_C1", "PLDE_C1", "POLE_C1", "PEAI_C1", "HYGL_C1", "TROR_C1", "TRCY_C1", "LARO_C1",
              "VERO_C2", "ARCA_C2", "PLDE_C2", "POLE_C2", "PEAI_C2", "HYGL_C2", "TROR_C2", "TRCY_C2", "LARO_C2",
              "VERO_C3", "ARCA_C3", "PLDE_C3", "POLE_C3", "PEAI_C3", "HYGL_C3", "TROR_C3", "TRCY_C3", "LARO_C3",
              "VERO_E1", "ARCA_E1", "PLDE_E1", "POLE_E1", "PEAI_E1", "HYGL_E1", "TROR_E1", "TRCY_E1", "LARO_E1",
              "VERO_E2", "ARCA_E2", "PLDE_E2", "POLE_E2", "PEAI_E2", "HYGL_E2", "TROR_E2", "TRCY_E2", "LARO_E2",
              "VERO_E3", "ARCA_E3", "PLDE_E3", "POLE_E3", "PEAI_E3", "HYGL_E3", "TROR_E3", "TRCY_E3", "LARO_E3",
              "VERO_T1", "ARCA_T1", "PLDE_T1", "POLE_T1", "PEAI_T1", "HYGL_T1", "TROR_T1", "TRCY_T1", "LARO_T1",
              "VERO_T2", "ARCA_T2", "PLDE_T2", "POLE_T2", "PEAI_T2", "HYGL_T2", "TROR_T2", "TRCY_T2", "LARO_T2")
set.seed(82)
RandomisedVector <- sample(Subplots, 72, replace = FALSE)
MatrixVector <- matrix(RandomisedVector, nrow = 8, ncol = 9)
#View(MatrixVector)
write.csv(MatrixVector, file = "Output/Site_layouts_2020/Site_8_B.csv")

Subplots <- c("VERO_C1", "ARCA_C1", "PLDE_C1", "POLE_C1", "PEAI_C1", "HYGL_C1", "TROR_C1", "TRCY_C1", "LARO_C1",
              "VERO_C2", "ARCA_C2", "PLDE_C2", "POLE_C2", "PEAI_C2", "HYGL_C2", "TROR_C2", "TRCY_C2", "LARO_C2",
              "VERO_C3", "ARCA_C3", "PLDE_C3", "POLE_C3", "PEAI_C3", "HYGL_C3", "TROR_C3", "TRCY_C3", "LARO_C3",
              "VERO_E1", "ARCA_E1", "PLDE_E1", "POLE_E1", "PEAI_E1", "HYGL_E1", "TROR_E1", "TRCY_E1", "LARO_E1",
              "VERO_E2", "ARCA_E2", "PLDE_E2", "POLE_E2", "PEAI_E2", "HYGL_E2", "TROR_E2", "TRCY_E2", "LARO_E2",
              "VERO_E3", "ARCA_E3", "PLDE_E3", "POLE_E3", "PEAI_E3", "HYGL_E3", "TROR_E3", "TRCY_E3", "LARO_E3",
              "VERO_T1", "ARCA_T1", "PLDE_T1", "POLE_T1", "PEAI_T1", "HYGL_T1", "TROR_T1", "TRCY_T1", "LARO_T1",
              "VERO_T2", "ARCA_T2", "PLDE_T2", "POLE_T2", "PEAI_T2", "HYGL_T2", "TROR_T2", "TRCY_T2", "LARO_T2")
set.seed(83)
RandomisedVector <- sample(Subplots, 72, replace = FALSE)
MatrixVector <- matrix(RandomisedVector, nrow = 8, ncol = 9)
#View(MatrixVector)
write.csv(MatrixVector, file = "Output/Site_layouts_2020/Site_8_C.csv")

########### Repeating above method for 2021!
##Adding 100 to each plot from last year
#sites #1-8, A = 1, B = 2, C = 3. e.g. site 1 control (C) = 13.
## e.g. site 4 'B' treatment = 42

#Only 8 species now, so 9 subplots per species plot. Want 5 of these to be thinned, 4 to be unthinned.
#Using same notation as last year: C for competition/unthinned, E for environment/thinned.
#Can assign trait plants later once we see what has germinated.
#No POLE

Subplots2021 <- c("VERO_C1", "ARCA_C1", "PLDE_C1", "PEAI_C1", "HYGL_C1", "TROR_C1", "TRCY_C1", "LARO_C1",
                  "VERO_C2", "ARCA_C2", "PLDE_C2", "PEAI_C2", "HYGL_C2", "TROR_C2", "TRCY_C2", "LARO_C2",
                  "VERO_C3", "ARCA_C3", "PLDE_C3", "PEAI_C3", "HYGL_C3", "TROR_C3", "TRCY_C3", "LARO_C3",
                  "VERO_C4", "ARCA_C4", "PLDE_C4", "PEAI_C4", "HYGL_C4", "TROR_C4", "TRCY_C4", "LARO_C4",
                  "VERO_E1", "ARCA_E1", "PLDE_E1", "PEAI_E1", "HYGL_E1", "TROR_E1", "TRCY_E1", "LARO_E1",
                  "VERO_E2", "ARCA_E2", "PLDE_E2", "PEAI_E2", "HYGL_E2", "TROR_E2", "TRCY_E2", "LARO_E2",
                  "VERO_E3", "ARCA_E3", "PLDE_E3", "PEAI_E3", "HYGL_E3", "TROR_E3", "TRCY_E3", "LARO_E3",
                  "VERO_E4", "ARCA_E4", "PLDE_E4", "PEAI_E4", "HYGL_E4", "TROR_E4", "TRCY_E4", "LARO_E4",
                  "VERO_E5", "ARCA_E5", "PLDE_E5", "PEAI_E5", "HYGL_E5", "TROR_E5", "TRCY_E5", "LARO_E5")

##Adding only reduced number of subplots, 6 total per species per plot with blanks otherwise

Subplotssix2021 <- c("VERO_C1", "ARCA_C1", "PLDE_C1", "PEAI_C1", "HYGL_C1", "TROR_C1", "TRCY_C1", "LARO_C1",
                  "VERO_C2", "ARCA_C2", "PLDE_C2", "PEAI_C2", "HYGL_C2", "TROR_C2", "TRCY_C2", "LARO_C2",
                  "VERO_C3", "ARCA_C3", "PLDE_C3", "PEAI_C3", "HYGL_C3", "TROR_C3", "TRCY_C3", "LARO_C3",
                  "VERO_E1", "ARCA_E1", "PLDE_E1", "PEAI_E1", "HYGL_E1", "TROR_E1", "TRCY_E1", "LARO_E1",
                  "VERO_E2", "ARCA_E2", "PLDE_E2", "PEAI_E2", "HYGL_E2", "TROR_E2", "TRCY_E2", "LARO_E2",
                  "VERO_E3", "ARCA_E3", "PLDE_E3", "PEAI_E3", "HYGL_E3", "TROR_E3", "TRCY_E3", "LARO_E3",
                  "Blank", "Blank", "Blank", "Blank", "Blank", "Blank", "Blank", "Blank",
                  "Blank", "Blank", "Blank", "Blank", "Blank", "Blank", "Blank", "Blank",
                  "Blank", "Blank", "Blank", "Blank", "Blank", "Blank", "Blank", "Blank")
#Using the above subplotssix for Sites 2 onwards

set.seed(111)
RandomisedVector <- sample(Subplots2021, 72, replace = FALSE)
MatrixVector <- matrix(RandomisedVector, nrow = 8, ncol = 9)
write.csv(MatrixVector, file = "Output/Site_layouts_2021/Site_1_A_2021.csv")
set.seed(112)
RandomisedVector <- sample(Subplots2021, 72, replace = FALSE)
MatrixVector <- matrix(RandomisedVector, nrow = 8, ncol = 9)
write.csv(MatrixVector, file = "Output/Site_layouts_2021/Site_1_B_2021.csv")
set.seed(113)
RandomisedVector <- sample(Subplots2021, 72, replace = FALSE)
MatrixVector <- matrix(RandomisedVector, nrow = 8, ncol = 9)
write.csv(MatrixVector, file = "Output/Site_layouts_2021/Site_1_C_2021.csv")
set.seed(121)
RandomisedVector <- sample(Subplotssix2021, 72, replace = FALSE)
MatrixVector <- matrix(RandomisedVector, nrow = 8, ncol = 9)
write.csv(MatrixVector, file = "Output/Site_layouts_2021/Site_2_A_2021.csv")
set.seed(122)
RandomisedVector <- sample(Subplotssix2021, 72, replace = FALSE)
MatrixVector <- matrix(RandomisedVector, nrow = 8, ncol = 9)
write.csv(MatrixVector, file = "Output/Site_layouts_2021/Site_2_B_2021.csv")
set.seed(123)
RandomisedVector <- sample(Subplotssix2021, 72, replace = FALSE)
MatrixVector <- matrix(RandomisedVector, nrow = 8, ncol = 9)
write.csv(MatrixVector, file = "Output/Site_layouts_2021/Site_2_C_2021.csv")
set.seed(131)
RandomisedVector <- sample(Subplotssix2021, 72, replace = FALSE)
MatrixVector <- matrix(RandomisedVector, nrow = 8, ncol = 9)
write.csv(MatrixVector, file = "Output/Site_layouts_2021/Site_3_A_2021.csv")
#Made a mistake for this one last year, ending up being 132 so now it's 232
set.seed(232)
RandomisedVector <- sample(Subplotssix2021, 72, replace = FALSE)
MatrixVector <- matrix(RandomisedVector, nrow = 8, ncol = 9)
write.csv(MatrixVector, file = "Output/Site_layouts_2021/Site_3_B_2021.csv")
set.seed(133)
RandomisedVector <- sample(Subplotssix2021, 72, replace = FALSE)
MatrixVector <- matrix(RandomisedVector, nrow = 8, ncol = 9)
write.csv(MatrixVector, file = "Output/Site_layouts_2021/Site_3_C_2021.csv")
set.seed(141)
RandomisedVector <- sample(Subplotssix2021, 72, replace = FALSE)
MatrixVector <- matrix(RandomisedVector, nrow = 8, ncol = 9)
write.csv(MatrixVector, file = "Output/Site_layouts_2021/Site_4_A_2021.csv")
set.seed(142)
RandomisedVector <- sample(Subplotssix2021, 72, replace = FALSE)
MatrixVector <- matrix(RandomisedVector, nrow = 8, ncol = 9)
write.csv(MatrixVector, file = "Output/Site_layouts_2021/Site_4_B_2021.csv")
set.seed(143)
RandomisedVector <- sample(Subplotssix2021, 72, replace = FALSE)
MatrixVector <- matrix(RandomisedVector, nrow = 8, ncol = 9)
write.csv(MatrixVector, file = "Output/Site_layouts_2021/Site_4_C_2021.csv")
set.seed(151)
RandomisedVector <- sample(Subplotssix2021, 72, replace = FALSE)
MatrixVector <- matrix(RandomisedVector, nrow = 8, ncol = 9)
write.csv(MatrixVector, file = "Output/Site_layouts_2021/Site_5_A_2021.csv")
set.seed(152)
RandomisedVector <- sample(Subplotssix2021, 72, replace = FALSE)
MatrixVector <- matrix(RandomisedVector, nrow = 8, ncol = 9)
write.csv(MatrixVector, file = "Output/Site_layouts_2021/Site_5_B_2021.csv")
set.seed(153)
RandomisedVector <- sample(Subplotssix2021, 72, replace = FALSE)
MatrixVector <- matrix(RandomisedVector, nrow = 8, ncol = 9)
write.csv(MatrixVector, file = "Output/Site_layouts_2021/Site_5_C_2021.csv")
set.seed(161)
RandomisedVector <- sample(Subplotssix2021, 72, replace = FALSE)
MatrixVector <- matrix(RandomisedVector, nrow = 8, ncol = 9)
write.csv(MatrixVector, file = "Output/Site_layouts_2021/Site_6_A_2021.csv")
set.seed(162)
RandomisedVector <- sample(Subplotssix2021, 72, replace = FALSE)
MatrixVector <- matrix(RandomisedVector, nrow = 8, ncol = 9)
write.csv(MatrixVector, file = "Output/Site_layouts_2021/Site_6_B_2021.csv")
set.seed(163)
RandomisedVector <- sample(Subplotssix2021, 72, replace = FALSE)
MatrixVector <- matrix(RandomisedVector, nrow = 8, ncol = 9)
write.csv(MatrixVector, file = "Output/Site_layouts_2021/Site_6_C_2021.csv")
set.seed(171)
RandomisedVector <- sample(Subplotssix2021, 72, replace = FALSE)
MatrixVector <- matrix(RandomisedVector, nrow = 8, ncol = 9)
write.csv(MatrixVector, file = "Output/Site_layouts_2021/Site_7_A_2021.csv")
set.seed(172)
RandomisedVector <- sample(Subplotssix2021, 72, replace = FALSE)
MatrixVector <- matrix(RandomisedVector, nrow = 8, ncol = 9)
write.csv(MatrixVector, file = "Output/Site_layouts_2021/Site_7_B_2021.csv")
set.seed(173)
RandomisedVector <- sample(Subplotssix2021, 72, replace = FALSE)
MatrixVector <- matrix(RandomisedVector, nrow = 8, ncol = 9)
write.csv(MatrixVector, file = "Output/Site_layouts_2021/Site_7_C_2021.csv")
set.seed(181)
RandomisedVector <- sample(Subplotssix2021, 72, replace = FALSE)
MatrixVector <- matrix(RandomisedVector, nrow = 8, ncol = 9)
write.csv(MatrixVector, file = "Output/Site_layouts_2021/Site_8_A_2021.csv")
set.seed(182)
RandomisedVector <- sample(Subplotssix2021, 72, replace = FALSE)
MatrixVector <- matrix(RandomisedVector, nrow = 8, ncol = 9)
write.csv(MatrixVector, file = "Output/Site_layouts_2021/Site_8_B_2021.csv")
set.seed(183)
RandomisedVector <- sample(Subplotssix2021, 72, replace = FALSE)
MatrixVector <- matrix(RandomisedVector, nrow = 8, ncol = 9)
write.csv(MatrixVector, file = "Output/Site_layouts_2021/Site_8_C_2021.csv")


  
  
  