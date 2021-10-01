
library(shiny)
library(tidyverse)
library(leaflet)
library(shinyBS)
library(shinycssloaders)
library(shinythemes)
library(highcharter)
library(leafgl)
library(tseries)
library(FactoMineR)
library(factoextra)

CountryLngLat2 = tibble()

# rsconnect::showLogs(streaming = TRUE)
cat(file=stderr(), "\n")

# Importing & Cleaning Database -------------------------------------------------------

dbQuakes = fst::read_fst("data/database_quakes.fst")
dbTectonic <- read_csv("data/tectonic.csv") %>% unique()
dbStudent <- read_csv("data/student.csv")

echantillon_temoin2 = tibble(x = rep(Inf,23232))
echantillon_temoin2$Latitude = runif(23232,min = -90,max=90)
echantillon_temoin2$Longitude = runif(23232,min = -180,max=180)

resStudent = t.test(dbStudent$`Distance quakes/tectonic_plates`,dbStudent$`Distance uniform_law/tectonic_plates`,alternative = "less")

dbQuakes <- dbQuakes %>% separate(Date, c("Month","Day","Year"), convert = TRUE,remove = FALSE)

#Removing wrong data
dbQuakes <-  dbQuakes %>% slice(-c(3379,7513,20651))



#Changing class of data 
dbQuakes$Date <- as.Date(dbQuakes$Date, "%m/%d/%Y")

dbVolcano = read_delim("data/database_volcano.csv", delim = ";")
dbVolcano$latitude = as.numeric(dbVolcano$latitude)
dbVolcano$elevation = as.numeric(dbVolcano$elevation)

CountryLngLat = read_csv("data/concap.csv") %>%
  select(CountryName,CapitalLatitude,CapitalLongitude)

dbEruptions = read_csv("data/database_eruptions.csv")

dbEruptions = dbEruptions %>% 
  mutate(StartDate = as.Date(paste0(start_year,"/",start_month,"/",start_day))) %>% 
  mutate(EndDate = as.Date(paste0(end_year,"/",end_month,"/",end_day))) %>% 
  mutate(eruptionsTime = as.numeric(difftime(EndDate, StartDate), units="days"))

nbEruptions = dbEruptions %>%
  group_by(volcano_number) %>%
  summarise("Number of Eruptions" = n(),
            "Time mean of Eruptions" = mean(eruptionsTime,na.rm = T))

dbVolcano = dbVolcano %>% left_join(nbEruptions, by = "volcano_number")

performance=function(age,perf) {
    db_kp=data.frame(age,perf)
    db_kp=na.omit(db_kp)
    perf = data.frame(db_kp %>% group_by(age) %>% summarise("Max VEI"=max(perf)))
    return(perf)
}
volcano_max_vei = performance(dbEruptions$volcano_name,dbEruptions$vei)

dbVolcano = left_join(dbVolcano,volcano_max_vei, by = c("volcano_name" = "age")) 

km2_per_country = read_csv2("data/database_km2_per_country.csv")
km2_per_country$`2018` = as.numeric(km2_per_country$`2018`)

dbTectonic$lat = as.numeric(dbTectonic$lat)
dbTectonic$lon = as.numeric(dbTectonic$lon)
dbTectonic$plate = str_to_upper(dbTectonic$plate)

world <- maps::map('world', fill=TRUE, col="transparent", plot=FALSE)
IDs <- sapply(strsplit(world$names, ":"), function(x) x[1])
world_sp <- maptools::map2SpatialPolygons(world, IDs=IDs,
                                          proj4string=sp::CRS("+proj=longlat +datum=WGS84"))
pointsSP <- sp::SpatialPoints(cbind(x = dbQuakes$Longitude, y= dbQuakes$Latitude),
                              proj4string=sp::CRS("+proj=longlat +datum=WGS84"))
indices <- sp::over(pointsSP, world_sp)
stateNames <- sapply(world_sp@polygons, function(x) x@ID)
dbQuakes$Country <- stateNames[indices]

pointsSP <- sp::SpatialPoints(cbind(x = dbEruptions$longitude, y= dbEruptions$latitude),
                              proj4string=sp::CRS("+proj=longlat +datum=WGS84"))
indices <- sp::over(pointsSP, world_sp)
stateNames <- sapply(world_sp@polygons, function(x) x@ID)
dbEruptions$Country <- stateNames[indices]


quake_country <- dbQuakes %>% filter(!is.na(Country))
eruptions_country <- dbEruptions %>% filter(!is.na(Country))
colorData <- quake_country[["Magnitude"]]
pal <- colorBin(palette = "YlOrRd", colorData, 7, pretty = TRUE)

map <- leaflet(options = leafletOptions(minZoom = 2, maxZoom = 18,preferCanvas = F )) %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addLayersControl(
    overlayGroups = "show tectonics plates",
    position = c("topleft"),
    options = layersControlOptions(collapsed = TRUE)
  ) %>% 
  addLegend("topright",
            pal=pal,
            values=colorData,
            title="Magnitude",
            layerId="colorLegend") %>%
  setMaxBounds(lng1 = -210, lng2 = 210, lat1 = -1200, lat2 = 1200)

#Creating tectonics plates boundaries
plates = unique(dbTectonic$plate)
for (plate1 in plates){
    
    data = dbTectonic[(dbTectonic$plate == plate1),] %>% unique() 
    BorneInf = 0
    
    for (i in 1:nrow(data)){
        
        if(i + 1 <= nrow(data)){
            
            if(abs(data$lon[i]-data$lon[i + 1])>300){
                
                table1 = data %>% slice((BorneInf + 1):i)
                BorneInf = i
                map = map %>% 
                    addPolylines(data = table1, lat = ~lat,lng =~lon,weight = 2,color = "grey", group = "show tectonics plates")
                
            }
            
        }
        
    }
    
    table1 = data %>%  slice((BorneInf + 1):nrow(data))
    map = map %>% addPolylines(data = table1, lat = ~lat,lng =~lon,weight = 2,color = "grey", group = "show tectonics plates")
    
}




#Data for Table
dbTable <- dbQuakes %>% select(Date,Latitude,Longitude,Type,Depth,Magnitude)

#Per Month per Year and per Type
sum_year <- dbQuakes %>% 
  select(Year,Month,Magnitude,Type,Country)  %>% 
  group_by(Year,Month,Type,Country)  %>% 
  summarize(count = mean(Magnitude)) %>% 
  mutate(`Date` = paste0(as.numeric(Year),"-",as.character(Month)))

sum_year$Date = sum_year$Date %>% zoo::as.yearmon()

sum_year_world <- dbQuakes %>% 
  select(Year,Month,Magnitude,Type)  %>% 
  group_by(Year,Month,Type)  %>% 
  summarize(count = mean(Magnitude)) %>% 
  mutate(`Date` = paste0(as.character(Year),"-",as.character(Month))) 

sum_year_world$Date = zoo::as.yearmon(sum_year_world$Date)

#Creation des moyennes glissantes
sum_year$meanAverage <- stats::filter(sum_year$count, rep(1,12), sides = 1)/12
sum_year_world$meanAverage <- stats::filter(sum_year_world$count, rep(1,12), sides = 1)/12



#Per Month and per Type
sum_month <- dbQuakes %>% 
  select(Month,Magnitude,Type)  %>% 
  group_by(Month,Type)  %>% 
  summarize(count = mean(Magnitude))

#Recodage de la variable Month
for (i in 1:12){
  sum_month$Month[sum_month$Month == i] <- month.name[i]
}

meanTimeTotal = mean(eruptions_country$eruptionsTime,na.rm = T)

sum_country_time_eruptions <- eruptions_country %>% 
  group_by(Country) %>% 
  summarise(Time_mean=mean(eruptionsTime,na.rm = T))

sum_country_earthquake <- quake_country %>%
  group_by(Country) %>%
  summarise(Earthquakes=n()) %>%
  arrange(desc(Earthquakes))

sum_country_eruptions <- eruptions_country %>%
  group_by(Country) %>%
  summarise(Eruptions=n()) %>% 
  arrange(desc(Eruptions))

mean_country_eruptions <- eruptions_country %>%
  group_by(Country) %>%
  summarise(VEI_mean=mean(vei,na.rm = TRUE)) %>% 
  arrange(desc(VEI_mean))

mean_country_earthquake <- quake_country %>%
  group_by(Country) %>%
  summarise(Magnitude_mean=mean(Magnitude,na.rm = TRUE)) %>% 
  arrange(desc(Magnitude_mean))

sum_country <- sum_country_earthquake %>% full_join(sum_country_eruptions, by = "Country")

sum_country = sum_country %>%
  left_join(km2_per_country, by = c("Country"="Country Name")) %>%
  left_join(mean_country_eruptions, by = "Country") %>% 
  left_join(mean_country_earthquake, by = "Country") %>% 
  left_join(sum_country_time_eruptions, by = "Country") %>% 
  mutate("Index Dangerosity" = (Earthquakes/`2018`*Magnitude_mean*100 + Eruptions/`2018`*VEI_mean*1000))

sum_country = sum_country %>% 
  mutate("Index Dangerosity Norm" = (`Index Dangerosity`-mean(`Index Dangerosity`, na.rm = T))/sd(`Index Dangerosity`, na.rm = T))


dbmclust = dbVolcano %>%
  select(`Number of Eruptions`,`elevation`,population_within_100_km)


res = hclust(dist(dbmclust), method="ward.D")

K=3
cl = cutree(res, K) 



acp = PCA(dbmclust, graph = F)
PlotACPInd = fviz_pca_ind(acp, habillage = as.factor(cl))
PlotACPVar = fviz_pca_var(acp)

dbmclust$clusters = cl %>% as.factor()

dbVolcano$tectonic_settings[dbVolcano$tectonic_settings == "Intraplate / Continental crust (>25 km)"] <- "Intraplate"
dbVolcano$tectonic_settings[dbVolcano$tectonic_settings == "Intraplate / Intermediate crust (15-25 km)"] <- "Intraplate"
dbVolcano$tectonic_settings[dbVolcano$tectonic_settings == "Intraplate / Oceanic crust (< 15 km)"] <- "Intraplate"

dbVolcano$tectonic_settings[dbVolcano$tectonic_settings == "Rift zone / Continental crust (>25 km)"] <- "Rift zone"
dbVolcano$tectonic_settings[dbVolcano$tectonic_settings == "Rift zone / Intermediate crust (15-25 km)"] <- "Rift zone"
dbVolcano$tectonic_settings[dbVolcano$tectonic_settings == "Rift zone / Oceanic crust (< 15 km)"] <- "Rift zone"

dbVolcano$tectonic_settings[dbVolcano$tectonic_settings == "Subduction zone / Continental crust (>25 km)"] <- "Subduction zone"
dbVolcano$tectonic_settings[dbVolcano$tectonic_settings == "Subduction zone / Crustal thickness unknown"] <- "Subduction zone"
dbVolcano$tectonic_settings[dbVolcano$tectonic_settings == "Subduction zone / Intermediate crust (15-25 km)"] <- "Subduction zone"
dbVolcano$tectonic_settings[dbVolcano$tectonic_settings == "Subduction zone / Oceanic crust (< 15 km)"] <- "Subduction zone"


dbVolcano$primary_volcano_type[dbVolcano$primary_volcano_type == "Caldera"] <- "Calderas"
dbVolcano$primary_volcano_type[dbVolcano$primary_volcano_type == "Caldera(s)"] <- "Calderas"

dbVolcano$primary_volcano_type[dbVolcano$primary_volcano_type == "Complex"] <- "Complexs"
dbVolcano$primary_volcano_type[dbVolcano$primary_volcano_type == "Complex(es)"] <- "Complexs"

dbVolcano$primary_volcano_type[dbVolcano$primary_volcano_type == "Lava cone"] <- "Lava cones"
dbVolcano$primary_volcano_type[dbVolcano$primary_volcano_type == "Lava cone(es)"] <- "Lava cones"
dbVolcano$primary_volcano_type[dbVolcano$primary_volcano_type == "Lava cone(s)"] <- "Lava cones"

dbVolcano$primary_volcano_type[dbVolcano$primary_volcano_type == "Lava dome"] <- "Lava domes"
dbVolcano$primary_volcano_type[dbVolcano$primary_volcano_type == "Lava dome(s)"] <- "Lava domes"

dbVolcano$primary_volcano_type[dbVolcano$primary_volcano_type == "Pyroclastic cone"] <- "Pyroclastic cones"
dbVolcano$primary_volcano_type[dbVolcano$primary_volcano_type == "Pyroclastic cone(s)"] <- "Pyroclastic cones"

dbVolcano$primary_volcano_type[dbVolcano$primary_volcano_type == "Shield"] <- "Shields"
dbVolcano$primary_volcano_type[dbVolcano$primary_volcano_type == "Shield(s)"] <- "Shields"


dbVolcano$primary_volcano_type[dbVolcano$primary_volcano_type == "Stratovolcano"] <- "Stratovolcanos"
dbVolcano$primary_volcano_type[dbVolcano$primary_volcano_type == "Stratovolcano(es)"] <- "Stratovolcanos"
dbVolcano$primary_volcano_type[dbVolcano$primary_volcano_type == "Stratovolcano?"] <- "Stratovolcanos"


dbVolcano$primary_volcano_type[dbVolcano$primary_volcano_type == "Tuff cones"] <- "Tuff cones"
dbVolcano$primary_volcano_type[dbVolcano$primary_volcano_type == "Tuff cone(s)"] <- "Tuff cones"

dbVolcano = dbVolcano %>%
  filter(tectonic_settings != "Unknown") %>% 
  filter(primary_volcano_type != "Crater rows" & 
           primary_volcano_type != "Compound" &
           primary_volcano_type != "Lava cones" &
           primary_volcano_type != "Maar(s)" &
           primary_volcano_type != "Pyroclastic shield" &
           primary_volcano_type != "Subglacial" &
           primary_volcano_type != "Tuff cones" & 
           primary_volcano_type != "Fissure vent(s)"
  ) %>% 
  filter(volcano_name != "Tacana")

CA = CA(table(dbVolcano$primary_volcano_type,dbVolcano$tectonic_settings),graph = F)
PlotCA = plot.CA(CA)
# UI ----------------------------------------------------------------------


ui <- navbarPage(
    title = "Tectonics plates analyse",
    
    theme = shinytheme("sandstone"),

    tabPanel("Home",
             icon = icon("globe", lib = "font-awesome"),
             
             includeCSS("www/styles.css"),
             tags$br(),
             fluidRow(
               column(2),
               column(8,
                 HTML('<h1 class="entry-title">Descriptive statistics </h2>'),
                 tags$div(id = "entry-title2"),
                 tags$br(),
                 tags$br(),
                 tags$br(),
                 selectInput("Country",
                   label = tags$h3("Pick a country :"),
                   choices = c("All Countries", unique(sum_country$Country))
                 ),
                ),
             ),
             
             fluidRow(
               column(2),
               column(8,
                      div(class = "wrap",
                          div(class = "one",
                              h2(textOutput("indic-1")),
                              tags$br(),
                              tags$p(
                                tags$b("Number of Earthquakes"),
                                tags$br(),
                                tags$span(style = "font-size:85%;",
                                          "from 1965 to 2016")
                                 )
                              ),
                          div(class = "two",
                              h2(textOutput("indic-2")),
                              tags$br(),
                              tags$p(
                                tags$b("Number of Eruptions"),
                                tags$br(),
                                tags$span(style = "font-size:85%;",
                                          "since Holocene")
                              )
                              ),
                          div(class = "three",
                              h2(textOutput("indic-3")),
                              tags$br(),
                              tags$p(
                                tags$b("Dangerosity Index"),
                                icon("question-circle", lib = "font-awesome",id = "question-1") %>% 
                                  tipify("Dangerosity Index = number of earthquakes per km²*magnitude mean + number of eruptions per km²*VEI mean ", placement="right", trigger = "hover"),
                                tags$br(),
                                tags$span(style = "font-size:85%;",
                                          "out of 10")
                              )
                             ),
                          div(class = "four",
                              h2(textOutput("indic-4")),
                              tags$br(),
                              tags$p(
                                tags$b("average time of eruptions "),
                                tags$br(),
                                tags$span(style = "font-size:85%;",
                                          "in days")
                              )
                          )
                          )
                     )
               
               
             ),
             tags$br(),
             tags$br(),
             fluidRow(
               column(2),
               column(4,
                      div(class = "wrap2",
                          fluidRow(
                            column(6,
                                   selectInput("serieT2-Choice",
                                    label = NULL,
                                    choices = c("Original data","Moving average on 12 month"))),
                            column(6
                                   
                                   )
                          ),
                          highchartOutput("serieT2") %>% withSpinner(color = "#a9e1d7"),
                          tags$br(),
                          actionButton("button", label = "Stationnary time series",icon("info-circle")))),
               column(4,
                      div(class = "wrap2",
                          tags$br(),
                          highchartOutput("serieT3") %>% withSpinner(color = "#a9e1d7"))),
               column(2)
               ),
             tags$br(),
             tags$br(),
             fluidRow(
               column(2),
               column(8,
                      div(class = "wrap3",
                          leafglOutput("mymap",height = 700),
                          verbatimTextOutput("nbObs")
                          )
                      )
               ),
             tags$br(),
             tags$br(),
             tags$br()
             

             ),
    
    tabPanel("Inferential test",
             fluidRow(
               column(2),
               column(8,
                      tags$br(),
                      HTML('<h1 class="entry-title">Clustering and statistical test </h2>'),
                      tags$div(id = "entry-title2"),
                      tags$br(),
                      tags$br(),
                      "On March 27, 1964, at 5:36 pm an earthquake of magnitude 9.2 occurred in Alaska.
                      It is the second-largest earthquake ever recorded, next to the M9.5 earthquake in Chile in 1960. This catastrophic event was a  great leap forward to the modern age of earthquake science. Most of what we know about earthquakes can be traced back to the geological research done after the great Alaskan earthquake.",
                      tags$br(),
                      "In this tab, we will try to understand some geological mechanisms with statistics. ",
                      tags$br(),
                      fluidRow(
                        column(6,
                               selectInput("dataset",
                                           label = tags$h3("Pick a dataset: "),
                                           choices = c("Volcanoes and Eruptions","Earthquakes")
                               ))
                      ),
                      tags$br(),
                      tags$br(),
                      tags$div(id = "entry-title3"),
                      conditionalPanel(
                        condition = "input.dataset == 'Volcanoes and Eruptions'",

                        
                        fluidRow(
                          column(2),
                          column(8,
                                 tags$br(),
                                 tags$br(),
                                 tags$h4("I. Principal Component Analysis & Hierarchical Clustering"),
                                 tags$br(),
                                 tags$br(),
                                 tags$p(class = "statis-text2",'In this section, we will try to summarize data with clustering and Principal Component Analysis(PCA). Then we will study the link beetwen primary volcano type and tectonics settings'),
                        tags$br(),
                        tags$p(class = "statis-text2","Principal component analysis is an unsupervised data analysis method that allows data visualization.")
                          )
                        ),
                        plotOutput("PlotACPVar",
                                   height = "500"
                                  ),
                        fluidRow(column(2),
                                 column(8,
                                        "As we can observe, elevation and population within 100 km seems to be correlated. We are now going to do ascending hierarchical clustering with the Ward.D method. ",
                                        tags$br(),
                                        tags$br(),
                                        selectInput("varCluster",
                                                    label = "Pick a variable : ",
                                                    choices = c("Number of Eruptions","elevation","population_within_100_km"))
                                        )),

                        highchartOutput("cluster"),
                        fluidRow(
                          column(2),
                          column(8)
                        ),
                        plotOutput("PlotACPInd",
                                   height = "500"),
                        tags$br(),
                        tags$br(),
                        fluidRow(
                          column(2),
                          column(8,tags$h4("II. Correspondence Analysis"))),
                        tags$br(),
                        fluidRow(
                          column(2),
                          column(8,
                                 "Here, we try to study the link between tectonics plate type (Subduction,Rift zone & Intraplate) with Volcanos type.
                                 The p-value of chi square independence test is very close to 0 (p-value = 2.828162e-39).
                                 That mean there is a significative link between the two variables. We can observe on this graph the link between the different modalities.")
                        ),
                        plotOutput("PlotCA",
                                   height = "500")
                      ),
                      conditionalPanel(
                        condition = "input.dataset == 'Earthquakes'",
                        tags$br(),
                        tags$p(class = "statis-text",
                               fluidRow(
                                 column(2),
                                 column(8,
                                        "The goal of this page is to study the relationship between worldwide earthquake distribution, and tectonic plate boundaries. For this analyse, we only select the earthquakes from our database. We create a control sample where latitude and longitude are randomly creates from an uniform law. We can observe scatters plots of our dataset and our control sample :  "),
                               )
                               ),
                        plotOutput("DistEarthquakes",
                                   height = "500"),
                        tags$br(),
                        tags$br(),
                        fluidRow(
                          column(2),
                          column(8,
                                 "Now we can calculate the distance between each earthquake and the nearest tectonic plate. ",
                                 plotOutput("StudentPlot"),
                                 "Let's compare mean with Student test. We do not need to verify normality (n1>30 & n2>30)",
                                 tags$br(),
                                 'Here is the p-value of our Student Test (alternative = "less") :',
                                 textOutput("Student"),
                                 tags$br(),
                                 "The test is significative and we reject (H_0), the earthquake sample has a mean significantly smaller than the test sample.",
                                 "That mean there is a link between eartquakes and tectonics plates.",
                                 tags$br(),
                                 tags$br(),
                                 tags$br(),
                                 tags$br(),
                                 tags$br(),
                        )
                      ))
                      
                      
               )
             
      ),
             icon = icon("search")
    ),
    
    tabPanel("Explore Data",
             dataTableOutput("speciesDataTable") %>% withSpinner(color = "#a9e1d7"),
             icon = icon("table")
    )
    
)  

server <- function(input, output, session) {
    
  
  output$speciesDataTable = renderDataTable(
    dbTable,
    options = list(
      scrollX = TRUE,
      pageLength = 10)
  )
    
  output$mymap <- renderLeafgl({
    map
  })  
  
  dbMap_reac = reactive({

    if(input$Country == "All Countries"){
      dbMap2 = quake_country
    } else {
      dbMap2 = quake_country %>% filter(Country %in% input$Country)
    }
    dbMap2
  })

    
    
  observe({
    
    if(input$Country == "All Countries"){
      
      leafletProxy("mymap") %>% 
        setView(lng = 30, lat = 30, zoom = 2)
      
    } else {
      
      CountryLngLat2 = CountryLngLat %>% filter(CountryName %in% input$Country) 
      
      leafletProxy("mymap") %>% 
        flyTo(lng = CountryLngLat2$CapitalLongitude[1], lat = CountryLngLat2$CapitalLatitude[1], zoom = 4)
    }
    
    colorData <- dbMap_reac()[["Magnitude"]]
    pal <- colorBin(palette = "YlOrRd", colorData, 7, pretty = TRUE)
    radius <- dbMap_reac()[["Magnitude"]]^6.2/ 2
    
    
    
    leafletProxy("mymap",data = dbMap_reac()) %>% 
      clearGroup(group = "data") %>%
      addCircles(~Longitude, ~Latitude,
                 radius=radius,
                 layerId=~Magnitude,
                 stroke=FALSE,
                 fillOpacity=0.4,
                 fillColor=~pal(colorData),
                 group = "data")
    

                 
  })

    DataInBounds <- reactive({
        if (is.null(input$mymap_bounds)){
          return(NULL)
        } else {
          bounds <- input$mymap_bounds
          latRng <- range(bounds$north, bounds$south)
          lngRng <- range(bounds$east, bounds$west)

          subset(dbMap_reac(),
                 Latitude >= latRng[1] & Latitude <= latRng[2] &
                   Longitude >= lngRng[1] & Longitude <= lngRng[2])
        }

    })

    observeEvent(input$mymap_shape_click,{
            click <- input$mymap_shape_click
            selected <- dbMap_reac()[dbMap_reac()$Latitude == click$lat & dbMap_reac()$Longitude == click$lng,]
            text <- as.character(tagList(
                paste0("Magnitude: ", selected$Magnitude), tags$br(),
                paste0("Year: ", selected$Year), tags$br(),
                paste0("Depth: ", selected$Depth), tags$br(),
                paste0("Type: ", selected$Type),tags$br(),
            ))
        leafletProxy("mymap") %>% clearPopups() %>%
            addPopups(click$lng, click$lat, text)

    })
    
    output$nbObs <- renderText({

        paste0("Number of earthquakes in map bounds: ",nrow(DataInBounds()))

    })
    
    # Datas 
    
    sum_year_reac = reactive({
      
      if(input$Country == "All Countries"){
        sum_year2 = sum_year_world
      } else {
        sum_year2 = sum_year %>% filter(Country %in% input$Country)
      }
      sum_year2
    })
    
    sum_country_reac = reactive({
      
      if(input$Country == "All Countries"){
        sum_country2 = sum_country
      } else {
        sum_country2 = sum_country %>% filter(Country %in% input$Country)
      }
      sum_country2
    })
    
    sum_month2 = reactive({
      
      if(input$Country == "All Countries"){
        sum_month2 = sum_month
      } else {
        sum_month2 = dbQuakes %>% 
          filter(Country %in% input$Country) %>% 
          select(Month,Magnitude,Type)  %>% 
          group_by(Month,Type)  %>% 
          summarize(count = mean(Magnitude))
      }
      for (i in 1:12){
        sum_month2$Month[sum_month2$Month == i] <- month.name[i]
      }
      sum_month2
    })
    
    
    # Indicators:
    
    output$`indic-1` <- renderText(
      
      if(input$Country == "All Countries")
      {
        sum(sum_country_reac()$Earthquakes, na.rm = T)
      }
      else
      {
        sum_country_reac()$Earthquakes %>% as.character()
        
      }
      
    )
    
    output$`indic-2` <- renderText(
      
      if(input$Country == "All Countries")
      {
        sum(sum_country_reac()$Eruptions,na.rm = T) %>% round(1)
      }
      else
      {
        sum_country_reac()$Eruptions %>% round(1) %>% as.character()
      }
      
    )
    
    output$`indic-3` <- renderText(
      
      if(input$Country == "All Countries")
      {
        "NA"
      }
      else
      {
        sum_country_reac()$`Index Dangerosity Norm` %>% round(1) %>% as.character()
      }
      
    )
    
    output$`indic-4` <- renderText(
      
      if(input$Country == "All Countries")
      {
        meanTimeTotal %>% round(1) %>% as.character()
      }
      else
      {
        sum_country_reac()$`Time_mean` %>% round(1) %>% as.character()
      }
      
    )
    
    # Highchart
    
    output$serieT2 <- renderHighchart({

                     
        if(input$`serieT2-Choice` == "Original data")
        {
          hchart(sum_year_reac(), "line",hcaes(x = Date, y = round(count,2), group = Type), lineWidth = 1) %>%
            hc_colors(c("#A6611A", "#DFC27D","#C9DCAF","#ECEAD3","red")) %>% 
            hc_yAxis(title = list(text = "Magnitude Average"))%>% 
            hc_title(text = paste("Magnitude mean in ",input$Country," since 1965"))

        } else 
          
        {
          hchart(sum_year_reac(), "line",hcaes(x = Date, y = round(meanAverage,2),group = Type), lineWidth = 1) %>%
            hc_colors(c("#A6611A", "#DFC27D","#C9DCAF","#ECEAD3")) %>% 
            hc_yAxis(title = list(text = "Moving Average"))%>%
            hc_title(text = paste("12-month moving average in ",input$Country," since 1965"))
        }
        
                     


      
    })
    
    output$serieT3 <- renderHighchart(
        hchart(sum_month2(), "line",hcaes(x = Month, y = round(count,2), group = Type), lineWidth = 1) %>% 
          hc_colors(c("#A6611A", "#DFC27D","#C9DCAF","#ECEAD3")) %>% 
          hc_yAxis(title = list(text = "Mean Magnitude"))%>% 
          hc_title(text = paste("Magnitude mean per month in ",input$Country," since 1965"))
        
      
    )
    
    stationary <- reactive({
      if(nrow(sum_year_reac())>30)
      {
        sum_year_ts = sum_year_reac()$count
        res = adf.test(sum_year_ts,alternative = "stationary")
        return(res$p.value)
      } 
      else 
      {
        return(NA)
      }

    })
    

    observe({
      if(is.na(stationary())){
        updateActionButton(
          session = getDefaultReactiveDomain(),
          "button",
          label = "Not enough data")  
      }
      else if(stationary()>0.05)
      {
        updateActionButton(
          session = getDefaultReactiveDomain(),
          "button",
          label = "Non stationary time series")  
      } 
      else 
      {
        updateActionButton(
          session = getDefaultReactiveDomain(),
          "button",
          label = "Stationary time series")  
      }
    })
    
    observeEvent(input$button, {

      
      
      showModal(modalDialog(
        title = "Augmented Dickey–Fuller Test",
        footer = NULL,
        easyClose = TRUE,
        "The augmented", tags$strong("Dickey–Fuller test (ADF) "),"tests the null hypothesis that a unit root is present in a time series sample.
        We test :",
        tags$br(),
        tags$br(),
        "(H0) : non stationnary time series ",tags$strong("against"), "  (H1) :  stationnary time series",
        tags$br(),
        tags$br(),
        "If a time serie is stationary, this mean there is no seasonal effects.",
        "Statistics calculated on the time series are also", tags$strong("consistents over time"),"like the mean or the variance.",
        tags$br(),
        tags$br(),
        "Here is the p-value of the Dickey–Fuller test : ", stationary(),
        tags$br(),
        if(is.na(stationary())){
          HTML("There is not enough data to do the Augmented Dickey–Fuller Test.")
        }
        else if(stationary()>0.05){
          HTML("The p-value is bigger than 0.05, the test is not significant, and we <strong>reject the alternative hypothesis</strong>. The time series is not stationary.")
        } else {
          HTML("The p-value is smaller then 0.05, the test is significant, and we <strong>do not reject the alternative hypothesis</strong>. The time series is significantly stationary.")
        }
      ))
    })
    
    
    output$PlotACPVar <- renderPlot(
      PlotACPVar
    )
    
    output$PlotACPInd <- renderPlot(
      PlotACPInd
    )
    
    output$PlotCA <- renderPlot(
      PlotCA 
    )
    
    output$DistEarthquakes <- renderPlot({
      par(mfrow=c(1,2))
      plot(x = dbQuakes$Longitude,
           y = dbQuakes$Latitude,
           xlab="Longitude",
           ylab="Latitude",
           main =expression("Distribution of earthquakes "*italic("(n=23,232)")))
      plot(x = echantillon_temoin2$Longitude,
           y = echantillon_temoin2$Latitude,
           xlab="Longitude",
           ylab="Latitude",
           main =expression("Distribution of control sample "*italic("(n=23,232)")))
    })
    
    
    output$cluster <- renderHighchart({
      
      x = dbmclust[[input$varCluster]]
      
      hcboxplot(x = x,
                var = dbmclust$clusters,
                outliers = FALSE,
                # col = dbmclust$clusters
                )%>%
         hc_chart(type = "column")
        # hc_colors(c("red","green","blue")) 
      
      
    })
    
    output$StudentPlot <- renderPlot({
      par(mfrow=c(1,2))
      boxplot(dbStudent$`Distance quakes/tectonic_plates`,
              main = "Distance between earthquakes and tectonics plates ",
              outline = FALSE)
      boxplot(dbStudent$`Distance uniform_law/tectonic_plates`,
              main = "Distance between control sample and tectonics plates ",
              outline = FALSE)
    })
    
    output$Student <- renderText({
      "p-value < 2.2e-16" 
      })
}

# Run the application 
shinyApp(ui = ui, server = server)

