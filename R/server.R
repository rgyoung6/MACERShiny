# Written by Rob Young at the University of Guelph in Ontario Canada, Sept, 2023
#********************************************Main program section***********************************************
################# Server Function #############################################

#server <- function(input, output, session) {
shinyAppServer <- function(input, output, session) {

  ############## Initialize variables ############################################

  #dada implement reactive values
  autoSeqDownloadDisplayString <- reactiveValues(data = NA)
  #Other Variables
  ASVFile <- reactiveValues(data = NA)
  provenanceDataFile <- reactiveValues(data = NA)
  output$ASVFileOut <- shiny::renderText({as.character("No data file selected")})
  output$provenanceDataFileOut <- shiny::renderText({as.character("No data file selected")})
  ASVFileTable<-reactiveValues(data = NA)
  provenanceDataFileTable<-reactiveValues(data = NA)

  # Get the path where all of the folders containing the fastq files are located
  volumes = shinyFiles::getVolumes()

  ################## AutoSeqDownload Submit Function #####################################
    #Connect this to the shinyChooseButton
    shinyFiles::shinyFileChoose(input, "autoSeqDownload", roots = volumes, session = session)

    # Get the file with the primer data for this analysis
    shiny::observeEvent(input$autoSeqDownload, {

      if(!is.null(input$autoSeqDownload)){
        tryCatch(
          expr = {
            autoSeqDownloadFile <- shinyFiles::parseFilePaths(volumes, input$autoSeqDownload)
            autoSeqDownloadDisplayString$data <- as.character(autoSeqDownloadFile$datapath)
            output$autoSeqDownloadDisplay <- shiny::renderText({as.character(autoSeqDownloadDisplayString$data)})
         },
         error = function(e){
           print("Error - Dada Location Button choose file cancelled")
           autoSeqDownloadLocation$data <- NA
         },
         warning = function(w){
           print("Warning - Dada Location Button choose file cancelled")
           autoSeqDownloadLocation$data <- NA
         }
       )
      }
    },ignoreInit = TRUE)

  shiny::observeEvent(input$dadaSubmit, {

   if (!is.na(autoSeqDownloadDisplayString$data) && is.character(autoSeqDownloadDisplayString$data) && length(autoSeqDownloadDisplayString$data) != 0){

     # Create local variables to call the function to ensure there are no errors

     input_file <- force(autoSeqDownloadDisplayString$data)
     if(is.na(input$searchStr) || is.null(input$searchStr) ){
       search_str<-NULL
     }else{
       search_str<-force(input$searchStr)
     }
     if(is.na(input$outputFile) || is.null(input$outputFile) ){
       output_file<-NULL
     }else{
       output_file<-force(input$outputFile)
     }
     NCBI_database <- force(input$ncbiInclude)
     BOLD_database <- force(input$boldInclude)
     seq_min <- force(input$autoSeqDownloadSeqMin)
     seq_max <- force(input$autoSeqDownloadSeqMax)

      tryCatch(
        expr = {
        #Run the Dada function here.

        shiny::showModal(shiny::modalDialog(
         title = "Automatic Sequence Download is underway.",
         "Processing, see the R console for details, please stand by...", footer=""

        ))

          #Run the function
          auto_seq_download(BOLD_database=BOLD_database,
                            NCBI_database=NCBI_database,
                            search_str= search_str,
                            input_file= input_file,
                            output_file=output_file,
                            seq_min=seq_min,
                            seq_max=seq_max)

         removeModal()
         shiny::showModal(shiny::modalDialog(
           title = "Automatic Sequence Download is complete",
           "Please see output files in the target directory."
         ))
       },
       error = function(e){
         removeModal()
         shiny::showModal(shiny::modalDialog(
           title = "Error",
           "AutoSeqDownlaod Genus Location Button choose file error. Please refer to the R
           consol for more information - 1"
         ))
         print("Error - AutoSeqDownlaod Genus Location Button choose file - 1")
         autoSeqDownloadLocation$data <- NA
       },
       warning = function(w){
         removeModal()
         shiny::showModal(shiny::modalDialog(
           title = "Warning",
           "AutoSeqDownlaod Genus Location Button choose file warning. Please refer to the R
           console for more information - 2"
         ))
         print("Warning - AutoSeqDownlaod Genus Location Button choose file - 2")
         autoSeqDownloadLocation$data <- NA
       })
    }else{
      shiny::showModal(shiny::modalDialog(
        title = "Missing Data",
        "Please select a proper Genus file and try submitting again!"
      ))
      print("Warning - AutoSeqDownlaod Genus Location Button choose file cancelled - 3")
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
           autoSeqDownloadLocation$data <- NA
         },
         warning = function(w){
           print("Warning - Dada Location Button choose file cancelled")
           autoSeqDownloadLocation$data <- NA
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

} # End of Server
