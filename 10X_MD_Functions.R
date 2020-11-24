#####################################################################
# Single Cell (10X) - QC Plot Aggregation Functions                 #
#####################################################################
# PROGRAM:                                                          #
#   10X_MD_Functions.R                                              #
# PROGRAMMER:                                                       #
#   Sarah Vang                                                      #
# DATE:                                                             #
#   11SEP2020                                                       #
# VERSION:                                                          #
#   R.3.5.0                                                         #
#-------------------------------------------------------------------#
# DESCRIPTION:                                                      #
#   This is a set of preparatory functions for the GEX / TCR 10X    #
#   Web Summary QC Compilation and assessment.                      #
#####################################################################

#-------------------------------------------------------------------#
# prepQC_Table:                                                     #
#   Performs some simple rearrangement for the QC metrics tables    #
#   compiled using 10X_QC_Compile.R.                                #
#                                                                   #
#      Params:                                                      #
#        - QC_Table  :                                              #
#            This is the read-in from the 10X_QC_Compile.R output.  #
#            It is a matrix (metric by sample) containing all of the#
#            QC Metrics for a single assay type.                    #
#-------------------------------------------------------------------#
prepQC_Table <- function(QC_Table){
  rownames(QC_Table) = QC_Table[,1]
  QC_Table = QC_Table[,c(-1)]
  return(QC_Table)
}

#-------------------------------------------------------------------#
# colorTable:                                                       #
#   This function takes in a QC dataframe from prepQC_Table and     #
#   produces a DT object (fancy RMD data table) which colors failed #
#   samples and bolds failed metrics. Thresholds for failure present#
#   in colorMani.                                                   #
#                                                                   #
#      Params:                                                      #
#        - dataFrame  :                                             #
#            This is data from the prepQC_Table function in a metric#
#            by sample data frame format.                           #
#        - colorMani  :                                             #
#            This is a data frame containing columns QC_Var,        #
#            Threshold, Direction. QC_Var must match the rows of    #
#            dataFrame and lists the metric under consideration.    #
#            Thresholds lists the numeric threshold of pass/fail    #
#            Direction (up/down) specifies whether failure occurs   #
#            above or below the threshold.                          #
#-------------------------------------------------------------------#
colorTable <- function(dataFrame,colorMani){
  if(!is.data.frame(dataFrame)){
    stop("ERROR - colorTable() - input data is not in a data frame format!")
  }
  if(!all(c("QC_Var","Threshold","Direction")%in%colnames(colorMani))){
    stop("ERROR - colorTable() - input data is missing required columns!")
  }
  if(!all(colorMani$QC_Var%in%rownames(dataFrame))|!all(rownames(dataFrame)%in%colorMani$QC_Var)){
    stop("ERROR - colorTable() - Color Manifest and input data frame are not compatible. Variable names must match!")
  }
  if(!all(tolower(colorMani$Direction)%in%c("up","down"))){
    stop("ERROR - colorTable() - Color Manifest is improperly formatted. {Direction} column must be populated by 'up' or 'down' only.")
  }
  
  ### Process the DataFrame
  revDF = apply(X=dataFrame,MARGIN=2,gsub,pattern="[%|,]",replacement="")
  flagViolation = matrix(0,nrow=nrow(dataFrame),ncol=ncol(dataFrame))
  
  for(i in 1:nrow(revDF)){
    currMet    = rownames(revDF)[i]
    dirDiag    = colorMani[which(colorMani$QC_Var==currMet),"Direction"]
    threshDiag = colorMani[which(colorMani$QC_Var==currMet),"Threshold"]
    if(dirDiag=="up"){
      flagViolation[i,] = ifelse(test = (revDF[i,]>threshDiag),yes = 1,no = 0)
    } else {
      flagViolation[i,] = ifelse(test = (revDF[i,]<threshDiag),yes = 1,no = 0)
    }
  }
  
  ### Set up new data-frame that includes indicators of "bad" QC metrics and the QC metric reports
  ### Will hide the indicators from view, but they will be present so that formatStyle() can use them
  ### to define attributes. 
  newDat = cbind(flagViolation,dataFrame)
  
  ### Grab samples with a violation
  badCols = which(colSums(flagViolation)>0)+ncol(dataFrame)
  
  
  ### Now create data table object and format
  DT::datatable(newDat,
                options = list(scrollX = TRUE,
                               fixedColumns = list(leftColumns=1,rightColumns=0),
                               columnDefs = list(list(visible=F,targets=1:ncol(flagViolation))))) %>%
    formatStyle(columns = badCols,
                backgroundColor = "#FF9999") %>%
    formatStyle(columns=1:ncol(dataFrame)+ncol(dataFrame),
                valueColumns = 1:ncol(dataFrame),
                fontWeight=styleEqual(c(1,0),c("bold","normal")),color = styleEqual(c(1,0),c("#0000FF","black"))) %>%
    formatStyle(columns = " ",
                backgroundColor = "#CCCCCC",fontWeight="bold")
}
