manduca <- function(num.interactions, corr.j) {

  #Defining the initial variables
  # Number of interactions == num.interactions
  # Correlation between traits == r
  # Number of generations
  generations = 1000
  #Population sizes
  size.plants=600
  size.animals=600
  #Mean of the traits
  media.z.inicial.i=0; media.op.inicial.i=0;
  media.z.inicial.j=0; media.op.inicial.j=0;
  #Standard deviation for the traits
  sd.z.inicial.i=1; sd.op.inicial.i=1
  sd.z.inicial.j=1; sd.op.inicial.j=1
  #Strength of competition
  VC=0.95
  #Alpha values for mutualism fitness component
  alpha.value.i=1; alpha.value.j=1;
  #Offspring error
  sd.random.segr.z=0.05; sd.random.segr.op=0.01


  #Keep values numeric
  num.interactions <- as.numeric(num.interactions);
  generations <- as.numeric(generations)

  #=====================================
  # Create dataframes to save data for branching plots
  # Species i output
  output.i<-as.data.frame(matrix("NA",0,8))
  colnames(output.i)<-c("generation","id","ecological.trait","op","sex","species","father","mother")
  output.i$generation <- as.numeric(output.i$generation)
  output.i$id <- as.numeric(output.i$id)
  output.i$ecological.trait <- as.numeric(output.i$ecological.trait)
  output.i$op <- as.numeric(output.i$op)
  output.i$sex <- as.factor(output.i$sex)
  output.i$species <- as.factor(output.i$species)
  output.i$father <- as.numeric(output.i$father)
  output.i$mother <- as.numeric(output.i$mother)
  #.....................................
  # Species j output
  output.j<-as.data.frame(matrix("NA",0,8))
  colnames(output.j)<-c("generation","id","ecological.trait","op","sex","species","father","mother")
  output.j$generation <- as.numeric(output.j$generation)
  output.j$id <- as.numeric(output.j$id)
  output.j$ecological.trait <- as.numeric(output.j$ecological.trait)
  output.j$op <- as.numeric(output.j$op)
  output.j$sex <- as.factor(output.j$sex)
  output.j$species <- as.factor(output.j$species)
  output.j$father <- as.numeric(output.j$father)
  output.j$mother <- as.numeric(output.j$mother)
  #.....................................
  # Defining the first step
  step <- 0
  #.....................................
  #Defining the factors of sex and species
  levels.sex <- factor(c('M','F'))
  levels.species <- factor(c("i", "j"))
  #.....................................
  p.generation=numeric(size.plants)
  p.id=numeric(size.plants)
  p.ecological.trait=numeric(size.plants)
  p.op=numeric(size.plants)
  p.sex=character(size.plants)
  p.species=character(size.plants)
  p.father=numeric(size.plants)
  p.mother=numeric(size.plants)
  #=====================================
  # Create individuals
  p.generation<-rep(step,size.plants)
  p.id<-seq(((step*1000)+1):((step*1000)+size.plants))
  p.sex <- as.vector(sample(levels.sex,size=size.plants,replace=TRUE))
  p.species <-rep("i",size.plants)
  p.father<-rep(0,size.plants)
  p.mother<-rep(0,size.plants)

  # Create traits
  # correlation matrix
  corr.i=1
  cor.mat.i<-matrix(c(1, corr.i, corr.i, 1), nrow=2)

  # covariance matrix
  traits.i = mvrnorm(n=size.plants, mu=c(media.z.inicial.i,media.op.inicial.i),
                     Sigma=cor2cov(cor.mat.i,sd=c(sd.z.inicial.i,sd.op.inicial.i)))

  ##Trait zA
  p.ecological.trait = traits.i[, 1]
  ##Trait oA
  p.op = traits.i[, 2]

  # Data.frame of species i
  i <-data.frame(generation=p.generation,id=p.id,ecological.trait=p.ecological.trait,
                 op=p.op,sex=p.sex,species=p.species,father=p.father,mother=p.mother)
  #.....................................
  q.generation=numeric(size.animals)
  q.id=numeric(size.plants)
  q.ecological.trait=numeric(size.animals)
  q.op=numeric(size.animals)
  q.sex=character(size.plants)
  q.species=character(size.animals)
  q.father=numeric(size.animals)
  q.mother=numeric(size.animals)
  #.....................................
  q.generation<-rep(step,size.animals)
  q.id<-seq(((step*1000)+1):((step*1000)+size.animals))
  q.sex<- as.vector(sample(levels.sex,size=size.animals,replace=TRUE))
  q.species<-rep("j",size.animals)
  q.father<-rep(0,size.animals)
  q.mother<-rep(0,size.animals)

  # Create traits of j species
  #correlation matrix
  cor.mat.j<-matrix(c(1, corr.j, corr.j, 1), nrow=2)

  #covariance matrix
  traits.j = mvrnorm(n=size.animals, mu=c(media.z.inicial.j,media.op.inicial.j),
                     Sigma=cor2cov(cor.mat.j,sd=c(sd.z.inicial.j,sd.op.inicial.j) ))

  ##Trait zA
  q.ecological.trait = traits.j[, 1]
  ##Trait oA
  q.op = traits.j[, 2]

  j <-data.frame(generation=q.generation,id=q.id,ecological.trait=q.ecological.trait,
                 op=q.op,sex=q.sex,species=q.species,father=q.father,mother=q.mother)
  #.....................................
  i <- transform(i,generation=as.numeric(generation),id=as.numeric(id), ecological.trait = as.numeric(ecological.trait), op = as.numeric(op))
  j <- transform(j, generation=as.numeric(generation),id=as.numeric(id),ecological.trait = as.numeric(ecological.trait), op = as.numeric(op))
  #.....................................
  i.tab<-as.matrix(i)
  write(t(i.tab),file=paste("output/output.i",num.interactions, corr.j,generations, ".Rtab",sep="_"),append=TRUE,sep=";",ncolumns=dim(i)[2]);
  j.tab<-as.matrix(j)
  write(t(j.tab),file=paste("output/output.j",num.interactions, corr.j, generations, ".Rtab",sep="_"),append=TRUE,sep=";",ncolumns=dim(j)[2]);
  #=====================================

  for (step in 1:generations){
    #=====================================
    # Intraspecific competition
    #calculate the difference between ecological traits for i species
    diff.z.intra.competition.i <- as.matrix(dist(i[,"ecological.trait"],upper=TRUE))
    #calculate alpha
    alpha.i <- exp(-VC*((diff.z.intra.competition.i^2)))
    mean.alpha.i <- (colSums(alpha.i[,])-1)/(length(alpha.i[1,])-1)
    #calculate the competition
    p.intra.competition <- 1-mean.alpha.i
    #join the competition to the data.frame
    plot.i <- cbind(i,p.intra.competition)
    #.....................................
    #the same for j species
    diff.z.intra.competition.j <- as.matrix(dist(j[,"ecological.trait"],upper=TRUE))
    alpha.j <- exp(-VC*((diff.z.intra.competition.j^2)))
    mean.alpha.j <- (colSums(alpha.j[,])-1)/(length(alpha.j[1,])-1)
    p.intra.competition <- 1-mean.alpha.j
    plot.j <- cbind(j,p.intra.competition)

    # Interespecific interactions
    # change the names to simplify the interactions
    plants <- plot.i
    animals <- plot.j
    num.plants <- length(plants$id)
    num.animals <- length(animals$id)
    # create a matrix of plants and animals
    plant.matrix <- matrix(plants$ecological.trait,num.plants,num.animals)
    animal.matrix <- t(matrix(animals$ecological.trait,num.animals,num.plants))
    # calculate the diference between traits of plants and animals
    dif.individuals <- abs(plant.matrix-animal.matrix)
    colnames(dif.individuals) <- animals$id
    rownames(dif.individuals) <- plants$id
    dif.individuals <- as.data.frame(dif.individuals)

    # Oreder for preferences of each individual (the more similar, the higher the preference)
    # Some plants can or not be chose
    preference.order <- sapply(dif.individuals,FUN=order)
    names.animals <- as.numeric(colnames(dif.individuals))
    #put ids to not mess with the real ids number
    standardized.ids <- seq(1,length(plants$id))
    index.animals <- seq(1,length(animals$id))
    index.animals <- cbind(animals$id,index.animals)
    ###
    # Intaracting animals
    animals.int <- sort((rep(names.animals,num.interactions)))
    # Which plants were chose by the animals
    plants.int <- as.vector(preference.order[seq(1:num.interactions),])
    ###

    interacting.plants <- plants[plants.int[],] #data.frame of interacting plants
    animals.int2 <- match(animals.int,index.animals[,1])
    interacting.animals <- animals[animals.int2[],] #data.frame of interacting animals

    # Data.frame of paired individuals
    paired.individuals <- cbind(interacting.plants,interacting.animals)
    colnames(paired.individuals)<-c("generation.plants","id.plants","ecological.trait.plants","op.plants","sex.plants","species.plants","father.plants","mother.plants","p.intra.competition.plants","generation.animals","id.animals","ecological.trait.animals","op.animals","sex.animals","species.animals","father.animals","mother.animals","p.intra.competition.animals")

    # Calculate the mutualistic component per step
    wbio.plants.step <- (exp(-alpha.value.i*(paired.individuals[,"ecological.trait.plants"]-paired.individuals[,"ecological.trait.animals"])^2))
    wbio.animals.step <- (exp(-alpha.value.j*(paired.individuals[,"ecological.trait.animals"]-paired.individuals[,"ecological.trait.plants"])^2))

    #sum interactions, when happen
    #wbio = Pmut
    wbio.plants <- aggregate(wbio.plants.step,by=list(id=paired.individuals$id.plants),FUN=sum)
    wbio.animals <- aggregate(wbio.animals.step,by=list(id=paired.individuals$id.animals),FUN=sum)
    wbio.plants <- 1+(wbio.plants)
    wbio.animals <- 1+(wbio.animals)

    # Remove the duplicates of plants
    mating.i <- paired.individuals[-which(duplicated(paired.individuals$id.plants)), ]
    mating.i <- subset(mating.i, select=c(1,2,3,4,5,6,7,8,9))

    # Remove the duplicates of animals
    if (num.interactions==1){
      mating.j <- paired.individuals
      mating.j <- subset(mating.j, select=c(10,11,12,13,14,15,16,17,18))
    }
    if (num.interactions>1){
      mating.j <- paired.individuals[-which(duplicated(paired.individuals$id.animals)), ]
      mating.j <- subset(mating.j, select=c(10,11,12,13,14,15,16,17,18))
    }

    # Order plants data.frame
    iii <- order(mating.i$id.plants)
    mating.i <- mating.i[iii,]

    #Adding the mutualism to the data.frame
    mating.i <- cbind(mating.i,wbio.plants$x)
    mating.j <- cbind(mating.j,wbio.animals$x)

    # Insert the plants who did not interacted
    interact <- mating.i$id.plants
    no.interact <- !plot.i$id%in%interact
    no.interact <- plot.i[which(no.interact==TRUE),]
    no.interact$wbio.plants <- 1
    colnames(no.interact) <- colnames(mating.i)
    mating.i <- rbind(mating.i,no.interact)
    #=====================================

    ## Divide males and females
    mating.i.males<- subset(mating.i, sex.plants=="M")
    mating.i.females<- subset(mating.i, sex.plants=="F")
    mating.j.males<- subset(mating.j, sex.animals=="M")
    mating.j.females<- subset(mating.j, sex.animals=="F")
    #.....................................
    # Available males
    available.males.i <- sum(mating.i.males$id.plants>0)
    available.males.j <- sum(mating.j.males$id.animals>0)
    #....................................
    # standardize the relative fitness
    #Competition fitness
    mating.i.males$p.intra.competition.plants <- mating.i.males$p.intra.competition.plants/max(mating.i.males$p.intra.competition.plants)
    mating.i.females$p.intra.competition.plants <- mating.i.females$p.intra.competition.plants/max(mating.i.females$p.intra.competition.plants)
    mating.j.males$p.intra.competition.animals <- mating.j.males$p.intra.competition.animals/max(mating.j.males$p.intra.competition.animals)
    mating.j.females$p.intra.competition.animals <- mating.j.females$p.intra.competition.animals/max(mating.j.females$p.intra.competition.animals)

    #Mutualism fitness
    mating.i.males$wbio.plants <- mating.i.males$wbio.plants/max(mating.i.males$wbio.plants)
    mating.i.females$wbio.plants <- mating.i.females$wbio.plants/max(mating.i.females$wbio.plants)
    mating.j.males$wbio.animals <- mating.j.males$wbio.animals/max(mating.j.males$wbio.animals)
    mating.j.females$wbio.animals <- mating.j.females$wbio.animals/max(mating.j.females$wbio.animals)
    #.....................................

    # Defining females mating probabilities
    #Competition*Mutualism
    p.mating.i.females<-mating.i.females$p.intra.competition.plants*mating.i.females$wbio.plants
    p.mating.i.females<-p.mating.i.females/sum(p.mating.i.females)
    #.....................................
    p.mating.j.females<-mating.j.females$p.intra.competition.animals*mating.j.females$wbio.animals
    p.mating.j.females<-p.mating.j.females/sum(p.mating.j.females)
    #.....................................
    p.mating.i.males<-mating.i.males$p.intra.competition.plants*mating.i.males$wbio.plants
    p.mating.i.males<-p.mating.i.males/sum(p.mating.i.males)
    #.....................................
    p.mating.j.males<-mating.j.males$p.intra.competition.animals*mating.j.males$wbio.animals
    p.mating.j.males<-p.mating.j.males/sum(p.mating.j.males)
    #.....................................
    # Generate i Offspring
    #.....................................
    o.generation <- numeric(size.plants);o.id<-numeric(size.plants);o.ecological.trait<-numeric(size.plants);o.op <- numeric(size.plants);o.sex <- numeric(size.plants);o.species <- numeric(size.plants);o.father<-numeric(size.plants);o.mother<-numeric(size.plants);
    #.....................................
    #########

    # Assortative mating choice for plants
    # Offspring without loop / i
    # Random error of offspring
    # Create error
    #matriz de correlacao
    cor.mat.i<-matrix(c(1, corr.i, corr.i, 1), nrow=2)
    #matriz de covariancia
    erros.i = mvrnorm(n=size.plants, mu=c(media.z.inicial.i,media.op.inicial.i),
                      Sigma=cor2cov(cor.mat.i,sd=c(sd.random.segr.z,sd.random.segr.op) ))
    ##Erro zA
    random.segr.z.i = erros.i[, 1]
    ##Erro oA
    random.segr.op.i = erros.i[, 2]
    #Mothers
    random.mother.i <- mating.i.females[sample(nrow(mating.i.females),size.plants,prob=p.mating.i.females,replace=TRUE),]
    # Mating choice by similar traits
    all.males.i.per.columns <- matrix( mating.i.males[,"ecological.trait.plants"],length(mating.i.males$id.plants),length(random.mother.i$id.plants))
    one.female.per.row <- t(matrix(t(random.mother.i[,"ecological.trait.plants"]),length(random.mother.i$id.plants),length(mating.i.males$id.plants)))
    choice <- as.data.frame(abs(one.female.per.row-all.males.i.per.columns))
    i.females.preference.vectors <- sapply(choice,FUN=order)
    choosen.males.vector <- i.females.preference.vectors[1,]
    standardized.male.positions <- seq(1,length(mating.i.males$id.plants))
    mating.i.males <- cbind(mating.i.males,standardized.male.positions)
    selected.males <- mating.i.males[as.vector(standardized.male.positions)%in%as.vector(choosen.males.vector), ]
    table.choosen.males <- as.data.frame(table(choosen.males.vector))
    selected.males.replicated <- selected.males[rep(rownames(selected.males),table.choosen.males$Freq),]
    choosen.males.vector <- as.data.frame(choosen.males.vector)
    #Fathers selected by females
    selected.fathers.plants <- selected.males.replicated[match(choosen.males.vector$choosen.males.vector,selected.males.replicated$standardized.male.positions),]
    #Sum of mothers and fathers traits
    sum.z<-random.mother.i$ecological.trait.plants+selected.fathers.plants$ecological.trait.plants
    sum.op<-random.mother.i$op.plants+selected.fathers.plants$op.plants

    #Offspring i
    o.generation <- rep(step,size.plants)
    o.id <- seq(((size.plants*step)+1),(size.plants*(step+1)))
    o.ecological.trait <- (sum.z/2)+random.segr.z.i
    o.op <- (sum.op/2)+random.segr.op.i
    o.sex <- as.vector(sample(levels.sex,size=size.plants,replace=TRUE))
    o.species <- rep("i",size.plants)
    o.father <- selected.fathers.plants$id.plants
    o.mother <- random.mother.i$id.plants
    offspring.i <- data.frame(generation=o.generation,id=o.id,ecological.trait=o.ecological.trait,op=o.op,sex=o.sex,species=o.species,father=o.father,mother=o.mother)

    #Fisherian runaway process for animals
    #Offspring without loop / j
    # Error of offspring
    #correlation matrix
    cor.mat.j<-matrix(c(1, corr.j, corr.j, 1), nrow=2)

    #covariance matrix
    erros.j = mvrnorm(n=size.animals, mu=c(media.z.inicial.j,media.op.inicial.j),
                      Sigma=cor2cov(cor.mat.j,sd=c(sd.random.segr.z,sd.random.segr.op) ))
    ##Erro zB
    random.segr.z.j = erros.j[, 1]
    ##Erro oB
    random.segr.op.j = erros.j[, 2]

    #Choose females according to probability of mating (Pmat)
    p.mating.mothers <- sample(mating.j.females$id.animals, size = size.animals, replace = T, prob = p.mating.j.females)
    p.mating.mothers <- sort(p.mating.mothers)

    #Choose who is going to be father
    # Create a matrix where the females are in the columns and the rows 5% of males
    matrix.choice<-matrix(0,nrow = round(table(mating.j.males$sex.animals)*0.05),
                          ncol = size.animals)

    ##Sample 5% of males to each female according to probability of mating (Pmat)
    matrix.choice <- matrix(data = (sample(x = mating.j.males$id.animals,
                                          size =dim(matrix.choice)[1]*dim(matrix.choice)[2],
                                          replace=TRUE, prob=p.mating.j.males)),
                            nrow = round(table(mating.j.males$sex.animals)*0.05),
                            ncol = size.animals)
    # Transpose matrix
    mothers.and.fathers <- as.data.frame(matrix.choice) %>% pivot_longer(cols = V1:V600,names_to = "mothers",
                                                  values_to = "fathers") # %>% filter(mothers == "V1")
    #Data frame with possible fathers
    possib.fathers <- right_join(mothers.and.fathers, mating.j.males %>% select(id.animals,op.animals), by = c("fathers" = "id.animals"))
    possib.fathers <- possib.fathers %>% tidyr::drop_na()

    #Select the fathers according to the maximum value of the sex trait
    selected.fathers.animals <- lapply(X = unique(possib.fathers$mothers), FUN = function(chose.mother){
      father <- possib.fathers %>% filter(mothers == chose.mother) %>%
                    select(fathers, op.animals) %>% filter(op.animals == max(op.animals)) %>%
                    select(fathers) %>% pull()
       return(data.frame(mother = chose.mother,father = father))
    }
    )

    #Selected parents of j species
    parents.j <- unique(ldply (selected.fathers.animals, data.frame))
    parents.j$mother <- p.mating.mothers

    # Data of mothers
    mothers.j <- as.data.frame(parents.j$mother)
    colnames(mothers.j) <- "id.animals"
    mothers.j <- left_join(mothers.j, mating.j.females, by = "id.animals")

    # Data of fathers
    fathers.j <- as.data.frame(parents.j$father)
    colnames(fathers.j) <- "id.animals"
    fathers.j <- left_join(fathers.j, mating.j.males, by = "id.animals")

    #Sum of mothers and fathers traits
    sum.z <- mothers.j$ecological.trait.animals+fathers.j$ecological.trait.animals
    sum.op <- mothers.j$op.animals+fathers.j$op.animals

    # Offspring j
    o.generation <- rep(step,size.animals)
    o.id <- seq(((size.animals*step)+1),(size.animals*(step+1)))
    o.ecological.trait <- (sum.z/2)+random.segr.z.j
    o.op <- (sum.op/2)+random.segr.op.j
    o.sex <- as.vector(sample(levels.sex,size=size.animals,replace=TRUE))
    o.species <- rep("j",size.animals)
    o.father <- fathers.j$id.animals
    o.mother <- mothers.j$id.animals
    offspring.j <- data.frame(generation=o.generation,id=o.id,ecological.trait=o.ecological.trait,
                             op=o.op,sex=o.sex,species=o.species,father=o.father,mother=o.mother)

    #########
    #.............................
    offspring.i->i;offspring.j->j
    names(offspring.i)->names.i
    names(offspring.j)->names.j
    offspring.i.tab<-as.matrix(offspring.i)
    write(t(offspring.i.tab),file=paste("output/output.i",num.interactions, corr.j, generations, ".Rtab",sep="_"),append=TRUE,sep=";",ncolumns=dim(offspring.i)[2])
    offspring.j.tab<-as.matrix(offspring.j)
    write(t(offspring.j.tab),file=paste("output/output.j",num.interactions, corr.j, generations, ".Rtab",sep="_"),append=TRUE,sep=";",ncolumns=dim(offspring.j)[2]);
    #.............................

    cat("nint=",num.interactions,"ksi=","current generation = ",step, "remaining =",generations-step,  "\n")
  }
}
