library(shinyWidgets)
library(flexdashboard)
library(shinyalert)
library(dplyr)

# Load data for tool
#########################
# Probability of exposure
dat <- readRDS("final_df_for_tool_exposure_vaxed_Bayes_2011_2020_18Nov2022.RDS")
rownames(dat) <- NULL

# Round estimates to 2 significant figures
dat$LL95CI <- signif(x = dat$LL95CI, digits = 2)
dat$UL95CI <- signif(x = dat$UL95CI, digits = 2)
dat$Median  <- signif(x = dat$Median, digits = 2)

dat <- rename(dat, Jurisdiction = State)

tib1 <- dat

######################
# Probability of death
dat <- readRDS("final_df_for_tool_death_vaxed_Bayes_2011_2020_18Nov2022.RDS")
rownames(dat) <- NULL

# Round estimates to 2 significant figures
dat$LL95CI <- signif(x = dat$LL95CI, digits = 2)
dat$UL95CI <- signif(x = dat$UL95CI, digits = 2)
dat$Median  <- signif(x = dat$Median, digits = 2)

dat <- rename(dat, Jurisdiction = State)

tib2 <- dat

# Sort the US states alphabetically so they look nice in the tool
us_states <- as.character(unique(dat$Jurisdiction))
us_states <- sort(us_states)

# Define UI for application that allows users to answer risk factor questions
shinyApp(
  ui = fluidPage(
    
    titlePanel(
      h1("Quantitative risk assessment for rabies", align = "center")
    ),
    
    div(img(src='shutterstock_1673657431-0_crop_smaller.png'), style="text-align: center;"), # Downloaded from Shutterstock. Contributor: Studio Ayutaka. https://www.shutterstock.com/image-vector/cute-woodland-animals-border-set-1673657431. License: Limited usage in print, advertising, and packaging. Unlimited web distribution.
    
    # Text explaining tool
    p("This Shiny app is a quantitative risk assessment tool for potential rabies virus exposures. 
     It applies to exposures that occurred in the US or Puerto Rico and does not
      apply to situations where the person was exposed to a bat but had no known contact with the bat."),
    
    p("This tool is intended for use by experienced public health professionals at the state and local level in the US only. 
      It is not intended for use by the general public.", 
      em("Please consult your local health department for advice on potential rabies virus exposures.")),
    
    
    h2("Purpose"),
    
    p("The aim of the tool is to promote judicious use of post-exposure prophylaxis (PEP).
      It estimates (i) the probability that an animal would test positive for rabies virus given that a person was exposed, Pr(rabid|exposure) 
      and (ii) the probability that a person would die from rabies given that they were exposed, Pr(death|exposure)."),
    
    h2("Instructions"),
    
    p("Select ONE response for each of the questions. If more than one response is relevant, choose the one with the higher risk 
      (i.e., if an individual was bitten on the head and leg, select head/neck)."),
    
    p("The animal categories are mutually exclusive. Before making a selection, read through all of the options
      as well as the descriptions in the methods section."),
    
    tabsetPanel(
      tabPanel("Pr(rabid|exposure)", 
               
               selectizeGroupUI(
                 id = "my-filters1",
                 inline = FALSE,
                 params = list(
                   
                   Animal = list(inputId =  "Animal", 
                                 title = "What is the species of animal responsible for the suspected rabies exposure?", 
                                 placeholder = 'select'),
                   Jurisdiction = list(inputId = "Jurisdiction", 
                                title = "In which jurisdiction did the exposure occur?", 
                                placeholder = 'select'),
                   Provoked = list(inputId = "Provoked", 
                                   title = "What circumstances led to the rabies exposure?", 
                                   placeholder = 'select'),
                   Healthy = list(inputId = "Healthy", 
                                  title = "What is the health status of the animal?", 
                                  placeholder = 'select'),
                   Vaccinated = list(inputId = "Vaccinated", 
                                     title = "Is the animal up to date on its rabies vaccination?", 
                                     placeholder = 'select')
                 )
               ),
               
               DT::dataTableOutput("table1"),
               
               h2("Interpretation"),
               
               h3("Pr(rabid|exposure)"),
               
               p("The estimated risk threshold is", strong("0.0004 (95% CI: 0.0002 – 0.0006)."),
                 "If Pr(rabid|exposure) < risk threshold, PEP may not be needed; however, if Pr(rabid|exposure) > risk threshold OR 
      Pr(rabid|exposure) = risk threshold, PEP may be indicated. Consider consulting your state health department
      or CDC rabies experts."),
               
      ),
      
      tabPanel("Pr(death|exposure)", 
               
               selectizeGroupUI(
                 id = "my-filters2",
                 inline = FALSE,
                 params = list(
                   Location_cat = list(inputId = "Location_cat", 
                                       title = "What is the location of the exposure?", 
                                       placeholder = 'select'),
                   Severity = list(inputId = "Severity", 
                                   title = "What is the severity of the exposure?", 
                                   placeholder = 'select'),
                   Animal = list(inputId =  "Animal", 
                                 title = "What is the species of animal responsible for the suspected rabies exposure?", 
                                 placeholder = 'select'),
                   Jurisdiction = list(inputId = "Jurisdiction", 
                                title = "In which jurisdiction did the exposure occur?", 
                                placeholder = 'select'),
                   Provoked = list(inputId = "Provoked", 
                                   title = "What circumstances led to the rabies exposure?", 
                                   placeholder = 'select'),
                   Healthy = list(inputId = "Healthy", 
                                  title = "What is the health status of the animal?", 
                                  placeholder = 'select'),
                   Vaccinated = list(inputId = "Vaccinated", 
                                     title = "Is the animal up to date on its rabies vaccination?", 
                                     placeholder = 'select')
                 )
               ),
               
               DT::dataTableOutput("table2"),
               
               h2("Interpretation"),
      
      h3("Pr(death|exposure)"),
      
      p("We can compare the Pr(death|exposure) to the probabilities of dying within one year from any cause.
      These probabilities are listed on the Social Security Administration website.
      We selected some rows from the 2019 period life table, which was used in the 2022 Trustees Report",
      a(href="https://www.ssa.gov/oact/STATS/table4c6.html",
        "https://www.ssa.gov"),":"),
      
      tags$ul(
        tags$li("10-year-old male: 0.000097"),
        tags$li("10-year-old female: 0.000095"),
        tags$li("20-year-old male: 0.0011"),
        tags$li("20-year-old female: 0.0004"),
        tags$li("30-year-old male: 0.0018"),
        tags$li("30-year-old female: 0.0008"),
        tags$li("40-year-old male: 0.0026"),
        tags$li("40-year-old female: 0.0014"),
        tags$li("50-year-old male: 0.0049"),
        tags$li("50-year-old female: 0.0030"),
        tags$li("60-year-old male: 0.011"),
        tags$li("60-year-old female: 0.007"),
        tags$li("70-year-old male: 0.022"),
        tags$li("70-year-old female: 0.015"),
        tags$li("80-year-old male: 0.056"),
        tags$li("80-year-old female: 0.041"),
        tags$li("90-year-old male: 0.16"),
        tags$li("90-year-old female: 0.13"),
        tags$li("100-year-old male: 0.35"),
        tags$li("100-year-old female: 0.30")),
      
      p("We can also compare the probabilities of dying generated by the tool to the lifetime probabilities of dying
      from other causes. The National Safety Council lists the lifetime odds of death for selected causes
      in the US in 2020 on their website", 
      a(href="https://injuryfacts.nsc.org/all-injuries/preventable-death-overview/odds-of-dying/",
        "https://injuryfacts.nsc.org"), ". We converted the odds to probabilities:"),
      
      tags$ul(
        tags$li("Cancer: 0.13"),
        tags$li("Opiod overdose: 0.015"),
        tags$li("Motor-vehicle crash: 0.0098"),
        tags$li("Drowning: 0.00098"),
        tags$li("Accidental gun discharge: 0.00013"),
        tags$li("Cataclysmic storm: 0.00003")),
      
      p("It is important to note that the Pr(death|exposure) for rabies is not a lifetime probability 
      like the ones listed above; the incubation period for rabies virus is typically 2-3 months.
      Thus, we would expect the probabilities from the National Safety Council to be even lower 
      if they were calculated over a shorter period, such as one year or one month."),
      
      )),
    
    # numericInput("value", label = "Type value to display risk level (green = low, yellow = moderate, red = high)", 
    #              min = 0, max = 1, value = 0.01, step = 0.000000001),
    # gaugeOutput("gauge"),
    
    h2("Methods"),
    
    h3("Animal categories"),
    
    p("Only mammals can transmit rabies virus. Animals were classified into risk groups for rabies based on 
      existing guidance as well as factors such as taxonomy, diet, and 
      geographical location. Here is the list of multi-species categories:"), 
    
    tags$ul(
      tags$li(strong("Small rodents, lagomorphs, and Eulipotyphla:"),
              "chipmunk, gopher, lagomorph, mole, mouse, muskrat, prairie dog, 
      rabbit, rat, rodent, shrew, squirrel, flying squirrel, vole"), 
      tags$li(strong("Rodents of unusual size:"),
              "beaver, groundhog, marmot, nutria, porcupine"), 
      tags$li(strong("Mesocarnivores: "),
              "coyote, badger, wolverine, otter, coati, lynx, ermine, 
      mink, fisher, marten, weasel, ringtail"),
      tags$li(strong("Equine: "),
              "horse, equine, donkey, mule,
      miniature donkey"),
      tags$li(strong("Hoofed animals: "),
              "deer, elk, moose, antelope, pronghorn, bison, Asian water buffalo, 
      reindeer, muskox, bighorn sheep, javelina"),
      tags$li(strong("Non-native wild: "),
              "African wild dog, binturong, kudu, serval, tapir, cheetah, rhinoceros, 
              tiger, kangaroo, zebra, 
            caracal, wallaby, leopard, giraffe, Siberian lynx, red panda, 
            llama, lion, chimpanzee, alpaca, jaguar, monkey,  
            macaque, gazelle, camel, camelid, kinkajou, yak, marmoset, 
            meerkat, capybara, hippopotamus, ocelot, oryx, eland, 
            lemur, hyrax, wallaroo, genet"),
      tags$li(strong("Large carnivores: "),
              "bear, cougar, sea lion, seal, wolf, mountain lion, panther, puma"),
      tags$li(strong("Pocket pets: "),
              "chinchilla, gerbil, guinea pig, hamster, hedgehog, sugar glider")),
    
    h3("NASPHV survey"),
    
    p("To estimate the risk threshold, a survey was conducted by the National 
      Association of State Public Health Veterinarians (NASPHV) among local and state public health officials.
      Respondents were asked whether PEP would be recommended given 24 scenarios which describe possible rabies 
      exposures in their state. Only the following information was provided:"),
    
    tags$ul(
      tags$li("Animal species/common name"), 
      tags$li("Whether the exposure was provoked"), 
      tags$li("Whether the animal is healthy"),
      tags$li("Whether the animal is up to date on its rabies vaccines (domestic animals only)")),
    
    p("We used logistic regression on the survey data to predict whether PEP would be recommended from the Pr(rabid|exposure)
      associated with each scenario."),
    
    h3("Data"),
    
    p("Parameter estimates for the location and severity of rabies exposures came from the literature.
      The positivity rates for animals as well as estimates regarding the circumstances
      of the exposure (provoked and health status) came from the national rabies surveillance 
      system in the US (10 years of data from 2011-2020). Jurisdiction-specific estimates were used for the positivity rate for 
      bats, cats, dogs, raccoons, skunks, and foxes; however, regional estimates were used
      for some jurisdictions when uncertainty was high."),
    
    p("Further details about the methods can be found in the accompanying paper (reference X)."),
    
    h2("References"),
    
    p("1. Babes V. Traité de la Rage. Paris: Librairie J.-B. Baillière et Fils; 1912. Chapter 6.
      Available from:", a(href="https://gallica.bnf.fr/ark:/12148/bpt6k5462676f.texteImage", 
                          "https://gallica.bnf.fr")),
    
    p("2. Baltazard M, Ghodssi M. Treatment of persons bitten by rabid wolves in Iran.
      Bull. Wld Hlth Org. 1954;10:797-803."),
    
    p("3. Ma X, Monroe BP, Wallace RM, et al. Rabies surveillance in the United States during 2019.
      JAVMA. 2021;258(11):1205-1220."),
    
    p("4. Shim E, Hampson K, Cleaveland S, Galvani AP. Evaluating the cost-effectiveness of rabies post-exposure
prophylaxis: A case study in Tanzania. Vaccine. 2009;27(51):7167-7172."),

   p("5. Vaidya SA, Manning SE, Dhankhar P, et al. Estimating the risk of rabies 
      transmission to humans in the U.S.: a Delphi analysis. BMC Public Health. 2010;10:278."),

  #p("6. Pieracci, E. et al. Vital Signs: Trends in Human Rabies Deaths and Exposures — United States, 
    #1938–2018. MMWR Morb Mortal Wkly Rep. 2009;68:524-528."),

h2("Image credit"),

p("Image downloaded from Shutterstock under a Standard License 
  (limited usage in print, advertising, and packaging. Unlimited web distribution.). 
  Contributor: Studio Ayutaka.", 
  a(href="https://www.shutterstock.com/image-vector/cute-woodland-animals-border-set-1673657431",
    "https://www.shutterstock.com")),

h2("Contact information"),

p("Maintainer: Kelly Charniga, kcharniga@cdc.gov")

  ),

server = function(input, output, session) {
  
  # Render the output for probability of exposure
  res_mod1 <- callModule(
    module = selectizeGroupServer,
    id = "my-filters1",
    data = tib1,
    vars = colnames(tib1)
  )
  
  output$table1 <- DT::renderDataTable({
    res_mod1()
  })
  
  # Render the output for probability of death table
  res_mod2 <- callModule(
    module = selectizeGroupServer,
    id = "my-filters2",
    data = tib2,
    vars = colnames(tib2)
  )
  
  output$table2 <- DT::renderDataTable({
    res_mod2()
  })
  
  # Pop-up warning not to proceed if animal is available for testing or quarantine
  shinyalert(
    title = "Warning",
    text = "Outcomes of public health evaluation or rabies virus testing of an animal supersede any risk calculations provided by this tool. Information provided here is intended to inform public health response when the health- or infection-status of a rabies-suspect animal cannot be confirmed by public health officials.",
    size = "s", 
    closeOnEsc = TRUE,
    closeOnClickOutside = FALSE,
    html = FALSE,
    type = "warning",
    showConfirmButton = TRUE,
    showCancelButton = FALSE,
    confirmButtonText = "OK",
    confirmButtonCol = "#64E03E",
    timer = 0,
    imageUrl = "",
    animation = TRUE
  )
  
  # output$gauge = renderGauge({
  #   gauge(input$value, 
  #         min = 0, 
  #         max = 1, 
  #         sectors = gaugeSectors(danger = c(0.1, 1), 
  #                                warning = c(0.01, 0.1),
  #                                success = c(0, 0.01)))
  # })
  
},

options = list(height = 500)
)
