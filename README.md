# Capstone I 

This was the first capstone project, of a two-part series for my chemical engineering senior design class. I learned a lot on this first project and as a result, I cleaned up the code massively when doing the second project. All of the code was in a single script, as well as the the code being not very modular. Since the capstone II project has the same concepts, but with a better code base. Please look at my [Capstone II project](https://github.com/wesleyZero/capstone_II/tree/main) to evaluate my coding abilities. For this project, I would like to simply write what I learned to do better in the second project. 

## What I learned (not) to do, and do in the second project 

1. I created a huge massive script, I (once again) vastly underestimated the scope of the project (how large it was going to become) so it became a massive mess. 

2. I put all of the files in the same directory (why???). I learned to put images in the /img folder and pdfs in the /pdf folder.

3. I did not use structures, and instead [used a TON of constants](https://github.com/wesleyZero/ChE_Capstone/blob/newBranch/Level3.m#L58-L370). At my internship at KARL STORZ I learned (from reading others code) that structures make things so much cleaner in MATLAB code. Why don't people just create classes and objects in MATLAB? I have no idea. The MATLAB-way seem to be to make a bunch of structures. 

4. Use [way too many lambda functions (anonymous functions)](https://github.com/wesleyZero/ChE_Capstone/blob/newBranch/Level3.m#L372-L444). I clearly like lambda functions and they are very nice for making some functions more readable. But the way I used them in this huge script wasn't making things easier to understand, it was making it harder. 


