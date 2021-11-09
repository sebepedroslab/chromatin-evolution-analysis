# functions

venn.two = function(list1,list2,catname1,catname2,main="Venn",
                    col1="green3",col2="magenta3",
                    eulerbool=T,print.mode=c("raw","percent")) {
  
  library(VennDiagram)
  library(grid)
  
  # compute set overlaps, intersections, etc
  list_uniq      = base::unique(c(list1, list2))
  list_intersect = base::intersect(list1, list2)
  list_diff_i    = base::setdiff(list1, list2)
  list_diff_j    = base::setdiff(list2, list1)
  
  # draw venn
  venn = draw.pairwise.venn(
    area1 = length(list1),
    area2 = length(list2),
    cross.area = length(list_intersect),
    category=c(catname1,catname2),
    fontfamily="Helvetica",
    cat.fontfamily = "Helvetica",
    col=c(col1,col2),ext.text = F,
    cat.col=c(col1,col2),
    ind = F, print.mode = print.mode,
    euler.d = eulerbool)
  grid::grid.newpage()
  grid::grid.draw(venn)
  grid::grid.text(main,x=0.5,y=0.95,
                  gp = gpar(fontsize = 8, fontface = "bold"))
  
  # output
  output = list(
    "catname_i"      = catname1,
    "catname_j"      = catname2,
    "title"          = main,
    "list_uniq"      = list_uniq,
    "list_intersect" = list_intersect,
    "list_diff_i"    = list_diff_i,
    "list_diff_j"    = list_diff_j,
    "venn" = venn
    )
  
  return(output)
  
}


venn.three = function(list1,list2,list3,
                      catname1,catname2,catname3,
                      main="Venn",
                      col1="green3",col2="magenta3",col3="orange",
                      eulerbool=T,print.mode=c("raw","percent")) {
  
  library(VennDiagram)
  library(grid)
  
  # compute set overlaps, intersections, etc
  intersect_n12  = base::intersect(list1, list2)
  intersect_n13  = base::intersect(list1, list3)
  intersect_n23  = base::intersect(list2, list3)
  intersect_n123 = base::intersect(intersect_n12, intersect_n13)
  
  # draw venn
  venn = draw.triple.venn(
    area1 = length(list1),
    area2 = length(list2),
    area3 = length(list3),
    n12 = length(intersect_n12),
    n13 = length(intersect_n13),
    n23 = length(intersect_n23),
    n123 = length(intersect_n123),
    category=c(catname1,catname2,catname3),
    fontfamily="Helvetica",
    cat.fontfamily = "Helvetica",
    col=c(col1,col2,col3),ext.text = F,
    cat.col=c(col1,col2,col3),
    ind = F, print.mode = print.mode,
    euler.d = eulerbool, scaled = T)
  grid::grid.newpage()
  grid::grid.draw(venn)
  grid::grid.text(main,x=0.5,y=0.95,
                  gp = gpar(fontsize = 8, fontface = "bold"))
  
  # output
  output = list(
    "venn" = venn,
    "catname_1"      = catname1,
    "catname_2"      = catname2,
    "catname_3"      = catname3,
    "intersect_12" = intersect_n12,
    "intersect_13" = intersect_n13,
    "intersect_23" = intersect_n23,
    "intersect_123" = intersect_n123
  )
  
  return(output)
  
}

