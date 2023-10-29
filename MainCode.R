#####
## Clustering starts at Line 42. You can directly start from Line 42 if you have the processed data (processedData.csv)
## set working direction
setwd("/Users/eric/Desktop/UCD/Courses/STA 260/Project")
source("STA260Rfunctions.R")

if(!file.exists("processedData.csv")){
    originalDataPath = "Africa_data.csv"
    augmentedDataPath = "Africa_data_7M.csv"
    
    ##### 
    ## Create neighborhood and sentence data
    sentencePath = "sentence.txt"
    if(!file.exists(sentencePath)){
        if(file.exists(augmentedDataPath)){
            sentence = formSentenceFromFile(inputPath = augmentedDataPath, outputPath = sentencePath)
        }else{
            largeData = CreateNeighborhood(inputPath = originalDataPath, outputPath = augmentedDataPath)
            sentence = formSentenceFromFile(inputPath = NULL, outputPath = sentencePath, data = largeData)
        }
    }
    
    ##### 
    ## IMPORTANT: train .py before running the following codes
    ## Run TSNE. 
    actorVectorsPath = "embeddingWord2Vec.vec"
    actorTSNEPath = "tsneEmbedding.csv"
    if(!file.exists(actorTSNEPath)){
        actorTSNE = vec2tsne(input_file = actorVectorsPath, output_file = actorTSNEPath, max_iter = 500)
    }else{
        actorTSNE = read.csv(actorTSNEPath)
    }
    
    ##### 
    ## Get country-scaled data with actor and duration infos
    data = matchVaribles(originalDataPath, augmentedDataPath, actorTSNE)
    somalia = getCountryData(data, getEventAfter = "01 Jan 2019")
    write.csv(somalia, "processedData.csv", row.names = FALSE)
}

##### 
## Run clustering
somalia = read.csv("processedData.csv") ## if you get error, please read and run all above line
runBaseline = FALSE
numClusters = 50

somaliaNormalized = myNormalize(somalia[, c(which(colnames(somalia) %in% c("latitude", "longitude")), grep("Embedding", colnames(somalia)), which(colnames(somalia) == "duration"))])
somaliaClustersObj = mySpecc(somaliaNormalized, numClusters, filePath = "somaliaClusters.RData", forceRun = FALSE)
somaliaClusters = factor(somaliaClustersObj@.Data)
contribution = calculateVariableContribution(somaliaNormalized, somaliaClusters)
plot(somalia$longitude, somalia$latitude, cex = 1, pch = 20, col = rainbow(max(as.numeric(somaliaClusters)))[somaliaClusters])
groupsNearMogadishu = getGroupsNear(somalia, Clusters = somaliaClusters, returnGroupSizeLargerThan = 1)
lapply(getDataFromGroups(somalia, somaliaClusters, groupsNearMogadishu), function(x){x[1:5, ]})
plotGroupsNear(somalia, somaliaClusters, groupsNearMogadishu)

if(runBaseline){
    ## Baseline 1
    somaliaClustersObjBaseline1 = mySpecc(somaliaNormalized[, 1:2], numClusters, filePath = "somaliaClustersBaseline1.RData", forceRun = FALSE)
    somaliaClustersBaseline1 = factor(somaliaClustersObjBaseline1@.Data)
    contributionBaseline1 = calculateVariableContribution(somaliaNormalized, somaliaClustersBaseline1)
    plot(somalia$longitude, somalia$latitude, cex = 1, pch = 20, col = rainbow(max(as.numeric(somaliaClustersBaseline1)))[somaliaClustersBaseline1], xlab = "", ylab = "")
    groupsNearMogadishu1 = getGroupsNear(somalia, Clusters = somaliaClustersBaseline1, returnGroupSizeLargerThan = 1)
    lapply(getDataFromGroups(somalia, somaliaClustersBaseline1, groupsNearMogadishu1), function(x){x[1:5, ]})
    plotGroupsNear(somalia, somaliaClustersBaseline1, groupsNearMogadishu1)
    
    ## Baseline 2
    somaliaClustersObjBaseline2 = mySpecc(somaliaNormalized[, c(1, 2, dim(somaliaNormalized)[2])], numClusters, filePath = "somaliaClustersBaseline2.RData", forceRun = FALSE)
    somaliaClustersBaseline2 = factor(somaliaClustersObjBaseline2@.Data)
    contributionBaseline2 = calculateVariableContribution(somaliaNormalized, somaliaClustersBaseline2)
    plot(somalia$longitude, somalia$latitude, cex = 1, pch = 20, col = rainbow(max(as.numeric(somaliaClustersBaseline2)))[somaliaClustersBaseline2], xlab = "", ylab = "")
    groupsNearMogadishu2 = getGroupsNear(somalia, Clusters = somaliaClustersBaseline2, returnGroupSizeLargerThan = 1)
    lapply(getDataFromGroups(somalia, somaliaClustersBaseline2, groupsNearMogadishu2), function(x){x[1:5, ]})
    plotGroupsNear(somalia, somaliaClustersBaseline2, groupsNearMogadishu2)
    
    ## Baseline 3
    somaliaClustersObjBaseline3 = mySpecc(somaliaNormalized[, -dim(somaliaNormalized)[2]], numClusters, filePath = "somaliaClustersBaseline3.RData", forceRun = FALSE)
    somaliaClustersBaseline3 = factor(somaliaClustersObjBaseline3@.Data)
    contributionBaseline3 = calculateVariableContribution(somaliaNormalized, somaliaClustersBaseline3)
    plot(somalia$longitude, somalia$latitude, cex = 1, pch = 20, col = rainbow(max(as.numeric(somaliaClustersBaseline3)))[somaliaClustersBaseline3], xlab = "", ylab = "")
    groupsNearMogadishu3 = getGroupsNear(somalia, Clusters = somaliaClustersBaseline3, returnGroupSizeLargerThan = 1)
    lapply(getDataFromGroups(somalia, somaliaClustersBaseline3, groupsNearMogadishu3), function(x){x[1:5, ]})
    plotGroupsNear(somalia, somaliaClustersBaseline3, groupsNearMogadishu3)
}

