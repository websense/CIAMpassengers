# CIAMpassengers
CIAM is a data driven approach for classifying long-term engagement of public transport riders at multiple temporal scales

This repository is the R code for "CIAM is a data driven approach for classifying long-term engagement of public transport riders at multiple temporal scales" by Rachel Cardell-Oliver and Doina Olaru.
The public transport passenger data for this project is private and so can not be shared here.  However, all the code used in the paper is provided so the study can be replicated with other datasets.

Abstract of the paper: Many human activities, including daily travel, show a mix of stable, intermittent and changing patterns in demand by individuals over time. However, the lack of continuous, long-term, passenger-linked data for public transport (PT) journeys means that we do not know how passenger ridership evolves in real-world networks. This paper proposes the CIAM model for classifying long-term passenger engagement with PT.  CIAM is a data-driven model combining year-on-year churn (C), monthly intensity (I), annual (A) and multi-year (M) engagement. Parameter search algorithms are used to ensure that the learned features are distinctive and robust. We evaluated CIAM using a 5-year dataset from a PT network with over 300 million journeys. CIAM identified distinct patterns of long-term ridership at multiple time scales. Although the total number of annual journeys was relatively stable over the five years, we found long-term differences between passenger subgroups.  Churn of passengers was a major factor in ridership with  only 55% of passengers retained from year to year. Patterns of annual engagement are often intermittent, so short-term snapshots of a few weeks are typically not good indicators for longer term engagement. Only 27% of high-frequency, full-fare riders still have the same level of engagement four years later compared with 55% who continue high-frequency engagement after only one year.

The codebase comprises:

GlobalConstants.R Load R libraries and data, define utility functions

ChurnYearonYear.R Calculate year-on-year churn of passengers in the dataset

IntensityMonthly.R Discover distinctive bands of monthly journey counts

AnnualEngagement.R Discover patterns of engagement (AE) with PT over a year

MultiYearEngagement.R Discover year to year transition probabilities between annual engagement classes

CIAMpassengerSegments.R Application using CIAM to select interesting subsets of passengers to support decision making for PT planning and operations
