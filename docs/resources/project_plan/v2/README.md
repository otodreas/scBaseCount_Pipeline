# Project Plan Submission v2

The file `Project_Plan.pdf` is the project plan.

I also included the file `project_specification_v2.pdf`, which I shared with Dag and Lotta, and which has the signature of my supervisor. The content of `project_specification_v2.pdf` and `Project_Plan.pdf` are more or less the same.

## Changes from submission v1

**The requested changes were**

1. The entire project hinges on: compare methods and see if they differ. There is no way to evaluate if the results are correct. I don't see any suggestion on how sides are right in case of conflict.
2. There is a circular evaluation problem, beause CyteType depends on the embedding structure. You are testing embedding &rarr; annotation pipeline, But evaluating using annotation outputs influenced by the same embedding.
3. Hypothesis is weakly falsifiable. You hypothesize that methods will differ in measurable ways. This is obvious: they operate completely differently. What are you testing here? That the methods are different or that they can measure the same data in the same way.
4. Overall, the hypothesis is:
  - Not biologically meaningful
  - Not risky
  - Not strongly falsifiable
  Even the “refutation” case (methods agree) is unlikely.
5. There is no formal statistical framework; there certainly is no biological validation layer (it's not mandatory, but it is added to the concerns of no validation).
6. The company needs to provide a statement that they agree that the outcome of this project will be published on GitHub.
7. The proposed project is too ambitious.

**The changes are addressed as follows**

1.
2.
3.
4.
5.
6. The file `project_specification_v2.pdf` (which was submitted in the first round) is signed by the CEO of Nygen and explicitly grants permission for the project to be published on GitHub (see §3.3)
7.