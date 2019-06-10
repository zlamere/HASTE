# Contributing to HATS
We welcome all kinds of input! Feel free to:
- Report a bug
- Discuss the state of the code
- Submit fixes
- Propose new features

## Writing Code for HATS
We're open to ANY contributions to improve the project, but expect *energetic* discussion to justify code not following these guidlines!
- HATS is written exclusively in modern Fortran (specifically Fortran 2008). We are not equipped to support a multi-language project at this time.
- Avoid the use of language extensions. Where extensions are used, specify conditional compilation (by a "C-style" preprocessor) for the code utilizing the language extension with alternative standard-compliant code as the default compilation path. Conditional flags for gFortran and Intel Fortran are currently in use in the code (GFORT and IFORT respectively).
- Strive to be self-documenting! Thoroughly comment your code (describe the procedure throughout the procedure) and use descriptive procedure/variable names.
- Avoid use of variables with the SAVE attribute (this includes module variables and some other data-sharing constructs). These can introduce challenges when adapting the code for massivley parallel architectures. Similarly, avoid use of recursive procedures.
- **NO deleted or depreciated language features.** This includes `GOTO`, `ASSIGN`, integer `FORMAT`, and `PAUSE` among others. Please use modern and standard-compliant code in these instances. Also avoid use of `COMMON` blocks and `DATA` statements (there are better ways to achieve the same behavior).
- **NO implicit typing... EVER.** All types **must** be explicitly declared. ALL modules and procedures **must** have **`Implicit None`** declared in their scope. Mixed-type arithmetic should be coded with explicit standard-compliant type-conversions.
- Order your indexes appropriately for efficent memory access. Fortran array indexes are in **column-major** order. If an array is going to be sliced or accessed sequentially frequently, the index along which it is going to be sliced should be the lowest dimension. i.e.:
  - Slicing:
    - `a(:,1)` accesses a slice of the array that is stored contiguously in memory (this is fast).
    - `a(1,:)` must construct the slice from elements that are **not** adjacent to one another in memory (this is slow).
  - Sequential access: The lowest dimension should be incremented the fastest so that the array is travesed contiguously in memory.
    ```
    Do i = 1,n
      Do j = 1,m
        x = a(j,i)
      End Do
    End Do
    ```
- Use consistent and standard-compliant style. (No tabs, 4-space indentation)

## Making Contributions
Pull requests are the best way to propose changes to the codebase. We actively welcome your pull requests:
- Fork the repository and create your branch from master.
- If you've added code that should be tested, add tests.
- Make sure your code lints.
- Issue a pull request.

## Any contributions you make will be under the GNU General Public License version 3 software license
In short, when you submit code changes, your submissions are understood to be under the same GPL-3.0 that covers the project. Feel free to contact the maintainers if that's a concern.

## Report bugs using Github's issues
Bugs may be submitted as issues on GitHub. A good issue tends to have:
- A quick summary and/or background
  - What you expected would happen
  - What actually happens
- Steps to reproduce
- Sample code if applicable, or links to relevant source lines
- Notes (possibly including why you think this might be happening, or stuff you tried that didn't work).
Be specific! You found the bug, you're submitting the bug: I hate to break it to you, but that makes you the world's expert on the bug at the moment.
