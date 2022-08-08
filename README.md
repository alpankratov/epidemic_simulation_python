Project Background
This project is a simple simulator of an epidemic inplemented in Python. In this model a person can have four states:
- Susceptible individuals can be infected with the disease.
- Recovered individuals are free of the disease and cannot pass it on to others. We assume that recovered individuals are not susceptible to the disease any more. We can also class individuals who are not susceptible to the disease (for example because they have built up immunity as a result of an infection with a related pathogen) as recovered, as they cannot be infected with the disease either.
- Infected individuals can spread the disease to susceptible individuals they are in contact with. They will either recover from the infection or die as a result of the infection.
- After individuals have died from the infection their state will not change.

Assume that the individuals live on a regular m * n grid, the first individual lives at (0; 0), the second at (0; 1), etc. with the (mn)-th individual living at (m - 1; n - 1).
The speed at which the epidemic spreads depends on how much contact there is between infected individuals and susceptible individuals.
This project assumes that infected individuals do not reduce their contact with others. This can be the case if symptoms are only mild for the majority of cases and/or if individuals become infectious before they develop symptoms (though the latter would alter the dynamic a little).