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
    main_directory <- ""
    # Split the current directory by "/"
    split_directory <- strsplit(current_directory, "/")
    for (i in seq_along(split_directory)) {
        if (split_directory[[i]] == "main") {
            main_directory <- paste(main_directory, split_directory[[i]], sep="/")
            break
        } else {
            main_directory <- paste(main_directory, split_directory[[i]], sep="/")
        }
    }
    
    return(main_directory)
}
