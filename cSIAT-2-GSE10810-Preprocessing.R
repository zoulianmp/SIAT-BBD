#-------------------------------------------------------------------------------
# R scripts for demonstration of the algorithms taught in the course SIAT-BBD
# Teacher: Fengfeng Zhou
# Email: FengfengZhou@gmail.com
# Web: http://healthinformaticslab.org/course/SIAT-BBD/
# Semester: October 2014
# Version 1.0.1
# Update: 2014-10-23
# Dataset: GSE10810
#   Development of Transcriptomic Biomarker Signature in Human Saliva to Detect
#   Lung Cancer
#   URL: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10810
# Restriction and disclaimer: 
#   These scripts are provided for the course project only. Please do not
#   re-distribute all these scripts, in part or completeness. For any other
#   purpose (including both commercial and non-profit purposes), please
#   contact the teacher for authorization and license.
#   
#   This section must be retained with all these course scripts.
#   
#   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#   SOFTWARE.
#-------------------------------------------------------------------------------

## Load packages
library(affy)   # Affymetrix pre-processing
library(limma)  # two-color pre-processing; differential expression
library(hgu133plus2cdf)


## import "phenotype" data, describing the experimental design
phenoData <- read.AnnotatedDataFrame("class-GSE10810.txt",sep = "")

## RMA normalization
celfiles <- "."
eset <- justRMA(phenoData=phenoData, celfile.path=celfiles, compress = TRUE)
egDataMatrix <- exprs(eset)
write.csv(egDataMatrix, file="class-matrix-GSE10810.csv")

## differential expression
combn <- factor(paste(pData(phenoData)[,1],pData(phenoData)[,2], sep = "_"))
design <- model.matrix(~combn) # describe model to be fit

fit <- lmFit(eset, design)  # fit each probeset to model
efit <- eBayes(fit)        # empirical Bayes adjustment
topTable(efit, number=15)      # table of top-ranked differentially expressed probesets
write.csv(efit, "class-result-GSE10810.csv")

## plot
par(mfrow=c(2,2))
plot(egDataMatrix[,1], egDataMatrix[,2])
abline(0,1)
plot(egDataMatrix[,2], egDataMatrix[,3])
abline(0,1)
plot(egDataMatrix[,3], egDataMatrix[,4])
abline(0,1)
plot(egDataMatrix[,4], egDataMatrix[,5])
abline(0,1)
par(mfrow=c(1,1))

