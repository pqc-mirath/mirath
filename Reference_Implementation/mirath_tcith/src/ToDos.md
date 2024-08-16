# To Do List

## Concerning EmulateMPC_μ

We need procedures performing the following

1. **MATRIX over FF** × **VECTOR over FFμ**
2. **VECTOR over FFμ** × **MATRIX over FF** 

Remark. We need to check if **Γ** belongs to FF or FFμ. If it is the former then **2.**, if not then we will need also

3. **VECTOR over FFμ** × **MATRIX over FFμ**

We can split the calculations of `(S + rndₛ)(C + rnd𝒸)` to avoid the "mixed" addition of matrices in FF and FFμ.
Then, we will need

4. **Matrix over FF** + **Matrix over FFμ**  
5. **Matrix over FF** × **Matrix over FF**


Additionally, we need to add the functions

6. `SplitCodeword()`
7. `Flatten()`

## Concerning EmulateParty_μ

7. **VECTOR over FF** × **Element in FFμ**

## Some comments/suggestions

For any function we add, we should include its respective tests.
Once all the tests successfully pass, then we should implement `EmulateMPC_μ()` and `EmulateParty_μ()`.

