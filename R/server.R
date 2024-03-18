# Written by Rob Young at the University of Guelph in Ontario Canada, Sept, 2023
#********************************************Main program section***********************************************
################# Server Function #############################################

#server <- function(input, output, session) {
shinyAppServer <- function(input, output, session) {

  ############## Initialize variables ############################################

  #dada implement reactive values
  primerFile <- reactiveValues(data = NA)
  #Other Variables
  ASVFile <- reactiveValues(data = NA)
  provenanceDataFile <- reactiveValues(data = NA)
  output$ASVFileOut <- shiny::renderText({as.character("No data file selected")})
  output$provenanceDataFileOut <- shiny::renderText({as.character("No data file selected")})
  ASVFileTable<-reactiveValues(data = NA)
  provenanceDataFileTable<-reactiveValues(data = NA)

  # Get the path where all of the folders containing the fastq files are located
  volumes = shinyFiles::getVolumes()

  ################## Dada Submit Function #####################################
    #Connect this to the shinyChooseButton
    shinyFiles::shinyFileChoose(input, "dadaDirectory", roots = volumes, session = session)

    # Get the file with the primer data for this analysis
    shiny::observeEvent(input$dadaDirectory, {

      if(!is.null(input$dadaDirectory)){
        tryCatch(
          expr = {

            dadaDirectory <- shinyFiles::parseFilePaths(volumes, input$dadaDirectory)

            dadaDirectoryDisplayString$data <- as.character(dadaDirectory$datapath)
            output$dadaDirectoryDisplay <- shiny::renderText({as.character(dadaDirectoryDisplayString$data)})

         },
         error = function(e){
           print("Error - Dada Location Button choose file cancelled")
           dadaLocation$data <- NA
         },
         warning = function(w){
           print("Warning - Dada Location Button choose file cancelled")
           dadaLocation$data <- NA
         }
       )
      }
    },ignoreInit = TRUE)

    #Connect this to the shinyChooseButton
    shinyFiles::shinyFileChoose(input, "primerFile", roots = volumes, session = session)

    # Get the file with the primer data for this analysis
    shiny::observeEvent(input$primerFile, {
      if(!is.null(input$primerFile)){
        tryCatch(
          expr = {

            primerFile <- shinyFiles::parseFilePaths(volumes, input$primerFile)
            primerFileDisplayString$data <- as.character(primerFile$datapath)
            output$primerFileDisplay <- shiny::renderText({as.character(primerFileDisplayString$data)})

          },
          error = function(e){
            print("Error - Primer Location Button choose file cancelled")
            dadaLocation$data <- NA
          },
          warning = function(w){
            print("Warning - Primer Location Button choose file cancelled")
            dadaLocation$data <- NA
          }
        )
      }
    },ignoreInit = TRUE)

  shiny::observeEvent(input$dadaSubmit, {

   if (!is.na(dadaDirectoryDisplayString$data) && is.character(dadaDirectoryDisplayString$data) && length(dadaDirectoryDisplayString$data) != 0 &&
       !is.na(primerFileDisplayString$data) && is.character(primerFileDisplayString$data) && length(primerFileDisplayString$data) != 0){

      # Create variables to call the dada_implement so that there are no conflicts
      # with the multithreading and the shiny app

      runFolderLoc <- force(dadaDirectoryDisplayString$data)
      primerFile <- force(primerFileDisplayString$data)

       if (force(input$uniOrbidirectional) == "Unidirectional"){
         unidirectional = TRUE
         bidirectional = FALSE
         fwdIdent <- ""
         revIdent <- ""
       }else if(force(input$uniOrbidirectional) == "Bidirectional"){
         unidirectional = FALSE
         bidirectional = TRUE
         fwdIdent <- force(input$fwdIdent)
         revIdent <- force(input$revIdent)
       }else{
         unidirectional = TRUE
         bidirectional = TRUE
         fwdIdent <- force(input$fwdIdent)
         revIdent <- force(input$revIdent)
       }
       printQualityPdf <- force(input$printQualityPdf)
       maxPrimeMis <- force(input$maxPrimeMis)
       fwdTrimLen <- force(input$fwdTrimLen)
       revTrimLen <- force(input$revTrimLen)
       maxEEVal <- force(input$maxEEVal)
       truncQValue <- force(input$truncQValue)
       truncLenValueF <- force(input$truncLenValueF)
       truncLenValueR <- force(input$truncLenValueR)
       error <- force(input$error)
       nbases <- force(input$nbases)
       maxMismatchValue <- force(input$maxMismatchValue)
       minOverlapValue <- force(input$minOverlapValue)
       trimOverhang <- force(input$trimOverhang)
       minFinalSeqLen <- force(input$minFinalSeqLen)

       tryCatch(
       expr = {
        #Run the Dada function here.

        shiny::showModal(shiny::modalDialog(
         title = "Dada analysis is underway.",
         "Processing, please stand by...", footer=""

        ))

        dada_implement(runFolderLoc = runFolderLoc,
                     primerFile = primerFile,
                     fwdIdent = fwdIdent,
                     revIdent = revIdent,
                     unidirectional = unidirectional,
                     bidirectional = bidirectional,
                     printQualityPdf = printQualityPdf,
                     maxPrimeMis = maxPrimeMis,
                     fwdTrimLen = fwdTrimLen,
                     revTrimLen = revTrimLen,
                     maxEEVal = maxEEVal,
                     truncQValue = truncQValue,
                     truncLenValueF = truncLenValueF,
                     truncLenValueR = truncLenValueR,
                     error = error,
                     nbases = nbases,
                     maxMismatchValue = maxMismatchValue,
                     minOverlapValue = minOverlapValue,
                     trimOverhang = trimOverhang,
                     minFinalSeqLen = minFinalSeqLen)
         removeModal()
         shiny::showModal(shiny::modalDialog(
           title = "Dada analysis is complete",
           "Please see output files in the target
           directory."
         ))
       },
       error = function(e){
         removeModal()
         shiny::showModal(shiny::modalDialog(
           title = "ERROR",
           "Dada Location Button choose file cancelled. Please refer to the R
           consol for more information - 1"
         ))
         print("Error - Dada Location Button choose file cancelled - 1")
         dadaLocation$data <- NA
       },
       warning = function(w){
         removeModal()
         shiny::showModal(shiny::modalDialog(
           title = "ERROR",
           "Dada Location Button choose file cancelled. Please refer to the R
           consol for more information - 2"
         ))
         print("Warning - Dada Location Button choose file cancelled - 2")
         dadaLocation$data <- NA
       })
    }else{
      shiny::showModal(shiny::modalDialog(
        title = "Missing Data",
        "Please select a primer file and try submitting again!"
      ))
      print("Warning - Dada Location Button choose file cancelled - 3")
    }
   },ignoreInit = TRUE)

  ################## Dada Combine Function ####################################
  #Get the location of the dada output files you would like to combine

    #Connect this to the shinyChooseButton
    shinyFiles::shinyFileChoose(input, "dadaCombineFile", roots = volumes, session = session)

    # Get the file with the primer data for this analysis
    shiny::observeEvent(input$dadaCombineFile, {

      if(!is.null(input$dadaCombineFile)){

        tryCatch(
          expr = {
            dadaCombineFile <- shinyFiles::parseFilePaths(volumes, input$dadaCombineFile)
            dadaCombineFileDisplayString$data <- as.character(dadaCombineFile$datapath)
            output$dadaCombineDisplay <- shiny::renderText({as.character(dadaCombineFileDisplayString$data)})
         },
         error = function(e){
           print("Error - Dada Location Button choose file cancelled")
           dadaLocation$data <- NA
         },
         warning = function(w){
           print("Warning - Dada Location Button choose file cancelled")
           dadaLocation$data <- NA
         }
       )
      }
    })

  #Running the data combine
  shiny::observeEvent(input$dadaCombine, {
    if (!is.na(dadaCombineFileDisplayString$data) && is.character(dadaCombineFileDisplayString$data) && length(dadaCombineFileDisplayString$data) != 0){
      # Create variables for the arguments to avoid conflicts between the multithreading
      # and the shiny
      fileLoc= force(dadaCombineFileDisplayString$data)
      minLen = force(input$dadaCombineMinLen)
      tryCatch(
        expr = {
          #Run the Dada function here.
          shiny::showModal(shiny::modalDialog(
            title = "Dada combine analysis results is underway.",
            "Processing, please stand by...", footer=""
          ))

        #Run the Dada combine function here.
          combine_dada_output(fileLoc = fileLoc, minLen = minLen)

        removeModal()
        shiny::showModal(shiny::modalDialog(
          title = "Dada combine analysis results is complete",
          "Please see output files in the target directory."
        ))
        },
        error = function(e){
          removeModal()
          shiny::showModal(shiny::modalDialog(
            title = "ERROR",
            "Please refer to the R consol for more information."
          ))
        },
        warning = function(w){
          removeModal()
          shiny::showModal(shiny::modalDialog(
            title = "ERROR",
            "Please refer to the R consol for more information."
          ))
        }
      )

    }else{
      shiny::showModal(shiny::modalDialog(
        title = "Missing Data",
        "Please select the file location where all of the DBTC dada output files are located that you wish to combine!"
      ))
    }
  },ignoreInit = TRUE)

  ################## Make BLAST DB Function ###################################

  #Connect this to the shinyChooseButton
  shinyFiles::shinyFileChoose(input, "makeBlastDBFileLoc", roots = volumes, session = session)

  # Get the fasta file you want to use to build your db
  shiny::observeEvent(input$makeBlastDBFileLoc, {
    tryCatch(
      expr = {

        makeBlastDBFileLoc <- shinyFiles::parseFilePaths(volumes, input$makeBlastDBFileLoc)
        makeBlastDBFileLocDisplayString$data <- as.character(makeBlastDBFileLoc$datapath)
        output$makeBlastDBFileLocDisplay <- shiny::renderText({as.character(makeBlastDBFileLocDisplayString$data)})

      },
      error = function(e){
        print("Error")
        makeBlastDBFileLoc$data <- NA
      },
      warning = function(w){
        print("Warning")
        makeBlastDBFileLoc$data <- NA
      }
    )
  },ignoreInit = TRUE)

  #Connect this to the shinyChooseButton
  shinyFiles::shinyFileChoose(input, "makeblastdbPath", roots = volumes, session = session)

  # Select where the makeblastdb program is on your computer
  shiny::observeEvent(input$makeblastdbPath, {
    tryCatch(
      expr = {

        makeblastdbPath <- shinyFiles::parseFilePaths(volumes, input$makeblastdbPath)
        makeblastdbPathDisplayString$data <- as.character(makeblastdbPath$datapath)
        output$makeblastdbPathDisplay <- shiny::renderText({as.character(makeblastdbPathDisplayString$data)})

      },
      error = function(e){
        print("Error")
        makeblastdbPath$data <- NA
      },
      warning = function(w){
        print("Warning")
        makeblastdbPath$data <- NA
      }
    )
  },ignoreInit = TRUE)

  #Connect this to the shinyChooseButton
  shinyFiles::shinyFileChoose(input, "makeBlastTaxaDBLoc", roots = volumes, session = session)

  # Point to the NCBI taxonomic data base
  shiny::observeEvent(input$makeBlastTaxaDBLoc, {
    tryCatch(
      expr = {

        makeBlastTaxaDBLoc <- shinyFiles::parseFilePaths(volumes, input$makeBlastTaxaDBLoc)
        makeBlastTaxaDBLocDisplayString$data <- as.character(makeBlastTaxaDBLoc$datapath)
        output$makeBlastTaxaDBLocDisplay <- shiny::renderText({as.character(makeBlastTaxaDBLocDisplayString$data)})

      },
      error = function(e){
        print("Error")
        makeBlastTaxaDBLoc$data <- NA
      },
      warning = function(w){
        print("Warning")
        makeBlastTaxaDBLoc$data <- NA
      }
    )
  },ignoreInit = TRUE)

  # Run the making BLAST db code.
  shiny::observeEvent(input$makeBlastDB, {
    if(is.na(makeblastdbPathDisplayString$data)){
      makeblastdbPathDisplayString$data <- "makeblastdb"
    }

    if (!is.na(makeBlastDBFileLocDisplayString$data) && is.character(makeBlastDBFileLocDisplayString$data) && length(makeBlastDBFileLocDisplayString$data) != 0 &&
      !is.na(makeBlastTaxaDBLocDisplayString$data) && is.character(makeBlastTaxaDBLocDisplayString$data) && length(makeBlastTaxaDBLocDisplayString$data)!= 0 ){

      # Create local variables to avoid conflicts with shiny and multithread
      fileLoc = force(makeBlastDBFileLocDisplayString$data)

      if(is.character(makeblastdbPathDisplayString$data) && length(makeblastdbPathDisplayString$data) != 0){

        makeblastdbPath = force(makeblastdbPathDisplayString$data)

      } else {

        makeblastdbPath = "makeblastdb"

      }

      taxaDBLoc = force(makeBlastTaxaDBLocDisplayString$data)
      dbName = force(input$dbName)
      minLen = force(input$makeBLASTDBMinLen)

      tryCatch(
        expr = {
          #Run the Dada function here.
          shiny::showModal(shiny::modalDialog(
            title = "Make BLAST database is underway.",
            "Processing, please stand by...", footer=""

          ))
          #Run the function
          make_BLAST_DB(fileLoc = fileLoc,
                        makeblastdbPath = makeblastdbPath,
                        taxaDBLoc = taxaDBLoc,
                        dbName = dbName,
                        minLen = minLen)
          removeModal()
          shiny::showModal(shiny::modalDialog(
            title = "Make BLAST database is complete",
            "Please see output files in the target
           directory."
          ))
        },
        error = function(e){
          removeModal()
          shiny::showModal(shiny::modalDialog(
            title = "ERROR",
            "There was an error in running the makeBLASTDB function. Make sure you have permissions for the target folder. Please see the R output for further details."
          ))
        },
        warning = function(w){
          removeModal()
          shiny::showModal(shiny::modalDialog(
            title = "ERROR",
            "There was an error in running the makeBLASTDB function. Make sure you have permissions for the target folder. Please see the R output for further details."
          ))
        })

    }else{
      shiny::showModal(shiny::modalDialog(
        title = "Missing Data",
        "Please fill in all of the necessary fields and submit again!"
      ))
    }
  },ignoreInit = TRUE)


  ################## BLAST sequences Function #################################

  #Connect this to the shinyChooseButton
  shinyFiles::shinyFileChoose(input, "BLASTDatabasePath", roots = volumes, session = session)

  # Point to the NCBI taxonomic data base
  shiny::observeEvent(input$BLASTDatabasePath, {
    tryCatch(
      expr = {
        BLASTDatabasePath <- shinyFiles::parseFilePaths(volumes, input$BLASTDatabasePath)
        BLASTDatabasePathDisplayString$data <- as.character(BLASTDatabasePath$datapath)
        output$BLASTDatabasePathDisplay <- shiny::renderText({as.character(BLASTDatabasePathDisplayString$data)})
      },
      error = function(e){
        print("Error")
        BLASTDatabasePath$data <- NA
      },
      warning = function(w){
        print("Warning")
        BLASTDatabasePath$data <- NA
      }
    )
  },ignoreInit = TRUE)

  #Connect this to the shinyChooseButton
  shinyFiles::shinyFileChoose(input, "blastnPath", roots = volumes, session = session)

  # Get the path where all of the folders containing the fastq files are located
  shiny::observeEvent(input$blastnPath, {
    tryCatch(
      expr = {
        blastnPath <- shinyFiles::parseFilePaths(volumes, input$blastnPath)
        blastnPathDisplayString$data <- as.character(blastnPath$datapath)
        output$blastnPathDisplay <- shiny::renderText({as.character(blastnPathDisplayString$data)})
      },
      error = function(e){
        print("Error")
        blastnPath$data <- NA
      },
      warning = function(w){
        print("Warning")
        blastnPath$data <- NA
      }
    )
  },ignoreInit = TRUE)

  #Connect this to the shinyChooseButton
  shinyFiles::shinyFileChoose(input, "querySeqPath", roots = volumes, session = session)

  # Point to the location of the fasta files you want to BLAST
  shiny::observeEvent(input$querySeqPath, {
    tryCatch(
      expr = {
        querySeqPath <- shinyFiles::parseFilePaths(volumes, input$querySeqPath)
        querySeqPathDisplayString$data <- as.character(querySeqPath$datapath)
        output$querySeqPathDisplay <- shiny::renderText({as.character(querySeqPathDisplayString$data)})
      },
      error = function(e){
        print("Error")
        querySeqPath$data <- NA
      },
      warning = function(w){
        print("Warning")
        querySeqPath$data <- NA
      }
    )
  },ignoreInit = TRUE)

  shiny::observeEvent(input$blastSequences, {

    if(is.na(blastnPathDisplayString$data) | is.null(blastnPathDisplayString$data)){
      blastnPathDisplayString$data <- "blastn"
    }

    if (!is.na(BLASTDatabasePathDisplayString$data) && is.character(BLASTDatabasePathDisplayString$data) && length(BLASTDatabasePathDisplayString$data) != 0 &&
        !is.na(blastnPathDisplayString$data) && is.character(blastnPathDisplayString$data) && length(blastnPathDisplayString$data) != 0 &&
        !is.na(querySeqPathDisplayString$data) && is.character(querySeqPathDisplayString$data) && length(querySeqPathDisplayString$data) != 0) {

      # Create local variables to avoid conflicts with shiny and multithread
      databasePath = force(BLASTDatabasePathDisplayString$data)
      blastnPath = force(blastnPathDisplayString$data)
      querySeqPath = force(querySeqPathDisplayString$data)
      minLen = force(input$BLASTminLen)
      BLASTResults = force(input$BLASTResults)
      numCores = force(input$blastSeqNumCores)

      tryCatch(
        expr = {
          #Run the Dada function here.

          shiny::showModal(shiny::modalDialog(
            title = "Sequence BLAST is underway.",
            "Processing, please stand by...", footer=""

          ))

          #Run the function
          seq_BLAST(databasePath = databasePath,
                    blastnPath = blastnPath,
                    querySeqPath = querySeqPath,
                    minLen = minLen,
                    BLASTResults = BLASTResults,
                    numCores = numCores)

          removeModal()
          shiny::showModal(shiny::modalDialog(
            title = "Sequence BLAST is complete",
            "Please see output files in the target directory."
          ))
        },
        error = function(e){
          removeModal()
          shiny::showModal(shiny::modalDialog(
            title = "ERROR",
            "Please refer to the R consol for more information."
          ))
        },
        warning = function(w){
          removeModal()
          shiny::showModal(shiny::modalDialog(
            title = "ERROR",
            "Please refer to the R consol for more information."
          ))
        }
      )
    }else{
      shiny::showModal(shiny::modalDialog(
        title = "Missing Data",
        "Please fill in all of the necessary fields and submit again!"
      ))
    }
  },ignoreInit = TRUE)

  ################## Taxon Assign Function ####################################

  #Connect this to the shinyChooseButton
  shinyFiles::shinyFileChoose(input, "taxaAssignFileLoc", roots = volumes, session = session)

  # Point to the BLAST output files to get taxonomic assignment
  shiny::observeEvent(input$taxaAssignFileLoc, {
    tryCatch(
      expr = {
        taxaAssignFileLoc <- shinyFiles::parseFilePaths(volumes, input$taxaAssignFileLoc)
        taxaAssignFileLocDisplayString$data <- as.character(taxaAssignFileLoc$datapath)
        output$taxaAssignFileLocDisplay <- shiny::renderText({as.character(taxaAssignFileLocDisplayString$data)})
      },
      error = function(e){
        print("Taxa Assign Error")
        taxaAssignFileLoc$data <- NA
      },
      warning = function(w){
        print("Taxa Assign Warning")
        taxaAssignFileLoc$data <- NA
      }
    )
  },ignoreInit = TRUE)

  #Connect this to the shinyChooseButton
  shinyFiles::shinyFileChoose(input, "taxaAssignDBLoc", roots = volumes, session = session)

  # Point to the NCBI taxonomic data base
  shiny::observeEvent(input$taxaAssignDBLoc, {
    tryCatch(
      expr = {

        taxaAssignDBLoc <- shinyFiles::parseFilePaths(volumes, input$taxaAssignDBLoc)
        taxaAssignDBLocDisplayString$data <- as.character(taxaAssignDBLoc$datapath)
        output$taxaAssignDBLocDisplay <- shiny::renderText({as.character(taxaAssignDBLocDisplayString$data)})

      },
      error = function(e){
        print("Error")
        taxaAssignDBLoc$data <- NA
      },
      warning = function(w){
        print("Warning")
        taxaAssignDBLoc$data <- NA
      }
    )
  },ignoreInit = TRUE)

  shiny::observeEvent(input$taxonAssign, {
    if (!is.na(taxaAssignFileLocDisplayString$data) && is.character(taxaAssignFileLocDisplayString$data) && length(taxaAssignFileLocDisplayString$data) != 0 && !is.na(taxaAssignDBLocDisplayString$data) && is.character(taxaAssignDBLocDisplayString$data) && length(taxaAssignDBLocDisplayString$data) != 0) {

      # Create local variables to avoid conflicts with shiny and multithread
       fileLoc = force(taxaAssignFileLocDisplayString$data)

       taxaDBLoc = force(taxaAssignDBLocDisplayString$data)
       numCores = force(input$taxaAssignNumCores)
       coverage = force(input$coverage)
       ident = force(input$ident)
       propThres = force(input$propThres)
       coverReportThresh = force(input$coverReportThresh)
       identReportThresh = force(input$identReportThresh)
       includeAllDada = force(input$includeAllDada)

       tryCatch(
         expr = {
           #Run the function here.
           shiny::showModal(shiny::modalDialog(
             title = "Taxonomic assingment is underway.",
             "See the R terminal for estimated time to completion.
             Processing, please stand by...", footer=""

           ))
           #Run the function
           taxon_assign(fileLoc = fileLoc,
                        taxaDBLoc = taxaDBLoc,
                        numCores = numCores,
                        coverage = coverage,
                        ident = ident,
                        propThres = propThres,
                        coverReportThresh = coverReportThresh,
                        identReportThresh = identReportThresh,
                        includeAllDada = includeAllDada)
            removeModal()

           shiny::showModal(shiny::modalDialog(
             title = "Taxonomic assingment is complete",
             "Please see output files in the target directory."
           ))
         },
         error = function(e){
           removeModal()
           shiny::showModal(shiny::modalDialog(
             title = "ERROR",
             "Please refer to the R consol for more information."
           ))
         },
         warning = function(w){
           removeModal()
           shiny::showModal(shiny::modalDialog(
             title = "ERROR",
             "Please refer to the R consol for more information."
           ))
         }
       )

    }else{
      shiny::showModal(shiny::modalDialog(
        title = "Missing Data",
        "Please fill in all of the necessary fields and submit again!"
      ))
    }
  })

  ################## Combine Taxa Assign Function #############################

  #Connect this to the shinyChooseButton
  shinyFiles::shinyFileChoose(input, "combineTaxaFileLoc", roots = volumes, session = session)

  # Point to the folder with the Taxon Assign records you would like to combine
  shiny::observeEvent(input$combineTaxaFileLoc, {
    tryCatch(
      expr = {

        combineTaxaFileLoc <- shinyFiles::parseFilePaths(volumes, input$combineTaxaFileLoc)
        combineTaxaFileLocDisplayString$data <- as.character(combineTaxaFileLoc$datapath)
        output$combineTaxaFileLocDisplay <- shiny::renderText({as.character(combineTaxaFileLocDisplayString$data)})

      },
      error = function(e){
        print("Error")
        combineTaxaFileLoc$data <- NA
      },
      warning = function(w){
        print("Warning")
        combineTaxaFileLoc$data <- NA
      }
    )
  },ignoreInit = TRUE)

  shiny::observeEvent(input$combineTaxa, {
    if (!is.na(combineTaxaFileLocDisplayString$data) && is.character(combineTaxaFileLocDisplayString$data) && length(combineTaxaFileLocDisplayString$data) != 0) {
      # Create local variables to avoid conflicts with shiny and multithread
      fileLoc = force(combineTaxaFileLocDisplayString$data)
      numCores = force(input$combineTaxaNumCores)

      tryCatch(
        expr = {
          #Run the function here.
          shiny::showModal(shiny::modalDialog(
            title = "Combining taxa assingment files is underway.",
            "Processing, please stand by...", footer=""

          ))
          #Run the function
          combine_assign_output(fileLoc = fileLoc,
                                numCores = numCores)
          removeModal()
          shiny::showModal(shiny::modalDialog(
            title = "Combining taxa assingment files is complete",
            "Please see output files in the target directory."
          ))
        },
        error = function(e){
          removeModal()
          shiny::showModal(shiny::modalDialog(
            title = "ERROR",
            "Please refer to the R consol for more information."
          ))
        },
        warning = function(w){
          removeModal()
          shiny::showModal(shiny::modalDialog(
            title = "ERROR",
            "Please refer to the R consol for more information."
          ))
        }
      )
    }else{
      shiny::showModal(shiny::modalDialog(
        title = "Missing Data",
        "Please fill in all of the necessary fields and submit again!"
      ))
    }
  })

  ################## Reduce Taxa Assign Function ##############################

  #Connect this to the shinyChooseButton
  shinyFiles::shinyFileChoose(input, "reduceTaxaFileLoc", roots = volumes, session = session)

  # Point to the folder with the Taxon Assign records you would like to combine
  shiny::observeEvent(input$reduceTaxaFileLoc, {
    tryCatch(
      expr = {

        reduceTaxaFileLoc <- shinyFiles::parseFilePaths(volumes, input$reduceTaxaFileLoc)
        reduceTaxaFileLocDisplayString$data <- as.character(reduceTaxaFileLoc$datapath)
        output$reduceTaxaFileLocDisplay <- shiny::renderText({as.character(reduceTaxaFileLocDisplayString$data)})

      },
      error = function(e){
        print("Error")
        reduceTaxaFileLocDisplayString$data <- NA
      },
      warning = function(w){
        print("Warning")
        reduceTaxaFileLocDisplayString$data <- NA
      }
    )
  },ignoreInit = TRUE)

  shiny::observeEvent(input$reduceTaxa, {
    if (!is.na(reduceTaxaFileLocDisplayString$data) && is.character(reduceTaxaFileLocDisplayString$data) && length(reduceTaxaFileLocDisplayString$data) != 0) {
      # Create local variables to avoid conflicts with shiny and multithread
      fileLoc = force(reduceTaxaFileLocDisplayString$data)
      numCores = force(input$reduceTaxaNumCores)

      tryCatch(
        expr = {
          #Run the function here.

          shiny::showModal(shiny::modalDialog(
            title = "Reduce ASV results with taxonomic assignment to unique taxa is underway.",
            "Processing, please stand by...", footer=""

          ))

          #Run the function
          reduce_taxa(fileLoc = fileLoc,
                      numCores = numCores)

          removeModal()
          shiny::showModal(shiny::modalDialog(
            title = "Reduce ASV results with taxonomic assignment to unique taxa is complete",
            "Please see output files in the target directory."
          ))
        },
        error = function(e){
          removeModal()
          shiny::showModal(shiny::modalDialog(
            title = "ERROR",
            "Please refer to the R consol for more information."
          ))
        },
        warning = function(w){
          removeModal()
          shiny::showModal(shiny::modalDialog(
            title = "ERROR",
            "Please refer to the R consol for more information."
          ))
        }
      )
    }else{
      shiny::showModal(shiny::modalDialog(
        title = "Missing Data",
        "Please fill in all of the necessary fields and submit again!"
      ))
    }
  })

  ################## Combine Reduce Taxa Assign Function ######################

  #Connect this to the shinyChooseButton
  shinyFiles::shinyFileChoose(input, "combineReducedTaxaFileLoc", roots = volumes, session = session)

  # Point to the folder with the Taxon Assign records you would like to combine
  shiny::observeEvent(input$combineReducedTaxaFileLoc, {
    tryCatch(
      expr = {

        combineReducedTaxaFileLoc <- shinyFiles::parseFilePaths(volumes, input$combineReducedTaxaFileLoc)
        combineReducedTaxaFileLocDisplayString$data <- as.character(combineReducedTaxaFileLoc$datapath)
        output$combineReducedTaxaFileLocDisplay <- shiny::renderText({as.character(combineReducedTaxaFileLocDisplayString$data)})

      },
      error = function(e){
        print("Error")
        combineReducedTaxaFileLocDisplayString$data <- NA
      },
      warning = function(w){
        print("Warning")
        combineReducedTaxaFileLocDisplayString$data <- NA
      }
    )
  },ignoreInit = TRUE)

  shiny::observeEvent(input$combineReduceTaxa, {

     if (!is.na(combineReducedTaxaFileLocDisplayString$data) && is.character(combineReducedTaxaFileLocDisplayString$data) && length(combineReducedTaxaFileLocDisplayString$data) != 0) {
      # Create local variables to avoid conflicts with shiny and multithread
      fileLoc = force(combineReducedTaxaFileLocDisplayString$data)
      presenceAbsence = force(input$presenceAbsence)

      tryCatch(
        expr = {
          #Run the function here.

          shiny::showModal(shiny::modalDialog(
            title = "Combine reduced taxonomic results for multiple markers for the same samples is underway.",
            "Processing, please stand by...", footer=""

          ))

          #Run the function
          combine_reduced_output(fileLoc = fileLoc, presenceAbsence = presenceAbsence)

          removeModal()
          shiny::showModal(shiny::modalDialog(
            title = "Combine reduced taxonomic results for multiple markers for the same samples is complete",
            "Please see output files in the target directory."
          ))
        },
        error = function(e){
          removeModal()
          shiny::showModal(shiny::modalDialog(
            title = "ERROR",
            "Please refer to the R consol for more information."
          ))
        },
        warning = function(w){
          removeModal()
          shiny::showModal(shiny::modalDialog(
            title = "ERROR",
            "Please refer to the R consol for more information."
          ))
        }
      )
    }else{
      shiny::showModal(shiny::modalDialog(
        title = "Missing Data",
        "Please fill in all of the necessary fields and submit again!"
      ))
    }
  })

} # End of Server
