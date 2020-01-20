# The GCG Pricing {#pricing}
> **This page is still in development and may be incomplete. Please excuse any inconveniences.**

## Pricing solvers
We currently have three different pricing solvers:

- MIP
- Knapsack
- Cliquer

All three of those can solve pricing problems heuristically or exact.\n

## Pricing Loop
The problem is then solved as follows:
- Loop through pricing subproblems in rounds
- Combine next pricing subproblem in queue with pricing solver
- Queue can be processed in parallel

## Best parameters for the pricing
We have around 2000 parameters for the pricing available, so it might be difficult
to choose the right ones for your class of problems. In the following, some hints
about good settings will be given.
