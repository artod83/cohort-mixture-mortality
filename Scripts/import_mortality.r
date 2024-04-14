m_list<-list.files("Life tables", pattern="*.csv")
m_years<-as.integer(gsub("[^0-9]+", "", m_list))
m<-lapply(m_list, function(x) {
  y<-read.csv(paste0("Life tables\\", x))
  if(ncol(y)==11){
    names(y)<-c("codreg",
                "region",
                "sex",
                "age",
                "lx",
                "dx",
                "qx",
                "Lx",
                "Px",
                "ex",
                "notes")
  } else {
    names(y)<-c("codreg",
                "region",
                "sex",
                "age",
                "lx",
                "dx",
                "qx",
                "Lx",
                "Px",
                "ex")    
  }
  y
  })
