Authors: Pietro Pollo, Shinichi Nakagawa & Michael Kasumovic

# Code Book

This code book describes all data fields in the `raw_data.csv` file (metadata). Each row in this file represents a effect size.


## Variables
- authors: authors of the empiric paper from which we extracted data
- year: year in which the empiric paper was published
- species: species used in experiments
- effect_id: identifier of effect sizes
- study_id: identifier of empirical studies used
- experiment_id: identifier of experiments
- measure_id: identifier of measure of male investment used
- species_id: identifier of species used
- f_trait: trait that varied in females in male mate choice experiments
- physical_contact: whether males could physically interact with females in experiments
- f_available: whether males had one (no-choice) or more females (multiple-choice) available in experiments
- f_quality_status: whether authors from the original paper presented evidence that female trait that varied in male mate choice experiments actually represented variation in benefits for males (i.e. female quality)
- m_trait: male trait evaluated
- m_trait_table1: which column from table 1 (see paper) the male trait evaluated belong to
- m_quality_status: whether authors from the original paper presented evidence that male trait evaluated represented variation in male capacity to acquire mates (i.e. male quality)
- m_description: male phenotype
- m_quality_cat: quality category of a male following information given by authors (if m_quality_status was verified)
- m_trait_cat: value trait category for certain traits (m_trait_table1 as column2)
- measure: type of male investment measured; whether pre-copula (latency and courtship), copula related or post-copulatory
- equation_type: how data was collected and the way effect sizes were calculated from it
- dependency: whether experimental measurements of male investment to high-quality females were dependent on male investment to low-quality females; important to calculate variance of effect sizes
- x1: mean male investment to high-quality females (equation_type = investment), proportion of male investment given high-quality females (equation_type = binary), or difference of male investment between high- and low-quality females (equation_type = difference); units vary
- x2: mean male investment to low-quality females (equation_type = investment), proportion of male investment given to low-quality females (equation_type = binary)
- sd1: standard deviation of x1
- sd2: standar deviation of x2
- n1: sample size of male investment to high-quality females
- n2: sample size of male investment to low-quality females
