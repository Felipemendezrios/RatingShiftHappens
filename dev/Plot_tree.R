library(data.tree)

# Your hierarchical data
hierarchy <- data.frame(
  indx = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14),
  level = c(1, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4),
  parent = c(0, 1, 1, 2, 2, 2, 3, 3, 3, 5, 5, 5, 7, 7),
  nS = c(2, 3, 3, 1, 3, 1, 2, 1, 1, 1, 1, 1, 1, 1)
)

hierarchy <- data.frame(
  indx = c(1,2,3,4,5,6,7,8,9,10,11,12),
  level = c(1,2,2,2,3,3,3,3,4,4,4,4),
  parent = c(0,1,1,1,2,2,3,3,5,5,7,7),
  nS = c(3,2,2,1,2,1,2,1,1,1,1,1)
)

hierarchy <- data.frame(
  indx = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16),
  level = c(1,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5),
  parent = c(0,1,1,2,2,2, 3,3, 4, 4,4,5,5,9,9,9),
  nS = c(2,3,2,3,2,1,1,1,3,1,1,1,1,1,1,1)
)

hierarchy <- data.frame(
  indx = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16),
  level = c(1,2,2,2,3,3,3,3,4,4,4,4,4,4,5,5),
  parent = c(0,1,1,1,2,2,4,4,5,5,6,6,8,8,12,12),
  nS = c(3,2,1,2,2,2,1,2,1,1,1,2,1,1,1,1)
)

# Function to find the path for a given node
find_path <- function(node, df) {
  if (node == 1) {
    return(as.character(node))
  }
  parent <- df[df$indx == node, "parent"]
  if(length(parent) > 0) {
    return(paste(find_path(parent, df), node, sep = "/"))
  } else {
    return(as.character(node))
  }
}

# Apply the function to all unique nodes to get their paths
paths <- sapply(hierarchy$indx, find_path, df = hierarchy)

# Create a data.tree structure from the paths
tree_data <- data.frame(pathString = paths)
tree <- as.Node(tree_data)

# Plot the tree structure
plot(tree)

SetGraphStyle(tree, rankdir = "TB")
SetNodeStyle(tree,style = "filled,rounded", shape = "box",
             fillcolor='red')

plot(tree)
