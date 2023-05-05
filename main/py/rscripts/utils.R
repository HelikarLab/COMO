# This file will hold various functions that will be useful throughout COMO

# Create a function that gets the current directory called "working_directory"
get_main_directory <-function() {
    #' This function will get the path of the "main" directory
    #' It will do this by getting the current working directory and then
    #'  splitting the path by "/". It will search through the items in the list
    #'  until one of the items equals "main". If the item is not main, it will append that item to a list
    #'  that holds the current path
    #'
    #' Returns: A string, "/home/username/main"

    current_directory <- getwd()
    main_directory <- c()
    
    # Split the current directory by "/", removing the first item
    split_directory <- unlist(strsplit(current_directory, "/"))[-1]
    
    
    for (item in split_directory) {
        # Append current directory to main_directory
        main_directory <- append(main_directory, item)
        if (item == "main") {
            break
        }
    }
    
    # Join main_directory by "/"
    main_directory <- paste(main_directory, collapse="/")
    main_directory <- paste0("/", main_directory)
    
    return(main_directory)
}
