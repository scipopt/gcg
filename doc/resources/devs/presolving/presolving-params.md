# Presolving Parameters {#presolving-params}
# Modifying Parameters of the GCG Presolving
There are different options to influence the presolving done to your problem. You can either
   * **Deactivate** presolving for your instance completely
   * **Change parameters** for the presolving

### Deactivate Presolving
The easiest way to not use presolving is to simply read in your instance and then, instead of doing an
`optimize`, perform `detect` and `optimize` instead. We will work on the unpresolved instance with its
respective decomposition. You will be able to read the following message in the log: 

```
GCG> optimize
there is an original decomposition and problem is not presolved yet -> disable presolving and start optimizing (rerun with presolve command before detect command for detecting in presolved problem)
```

### Change Parameters of the Presolving
You can also modify other parameters of the presolving. 
You can set them by entering `set` and then `presolving`. 
In the next section, we give a list of all parameters that can then be changed.
\n
A list of all presolving parameters can be seen (and searched) when inserting `set presolving` into the search box on the page of @ref interactive-menu.