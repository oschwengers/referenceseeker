# Contribution Guidelines
We welcome contributions!

We use the [GitHub Flow](https://guides.github.com/introduction/flow/) development style. Please set up a pull request against our master branch with any changes you want us to consider.

Additionally, it would be superb if the following conventions were followed:

## Check your coding style
- Make sure your contributions and changes follow [Python PEP 8](https://www.python.org/dev/peps/pep-0008/).
- Do not commit commented-out code or files that are no longer needed. Remove the code or the files unless there is a good reason to keep it.

## Use unit tests
To streamline the development process and to keep source codes in good shape, please:
- run the existing [PyTest](https://docs.pytest.org/en/latest/) unit tests before commiting new code:
```
$ cd <path-to-repository>
$ python3 -m pytest
```
- add tests to new code

## Make sure that your branch contains clean commits
- Follow the common sense guidelines for writing good commit messages (see below).
- Make separate commits for separate changes. If you cannot describe what the commit does in one sentence, it is probably a mix of changes and should be separated into several commits.
- Do not merge `master` into your branch. Use `git rebase` if you need to resolve merge conflicts or include the latest changes.

### Guidelines for good commit messages
1. Separate subject from body with a blank line
2. Use the imperative mood in the subject line ("Fix", "Add", "Change" instead of "Fixed", "Added", "Changed")
3. Limit the subject line to 50 characters
4. Reference an issue in the beginning of a subject line
5. Do not end the subject line with a period
6. Use the body to explain what and why vs. how

Bad:
<pre>
Some changes fixing this and that...
</pre>

Good:
<pre>
#123 Fix broken resistance model link in ABR reports

As CARD changed some of their website structure we need to update link templates for antibiotic resistance (ABR) models.
</pre>

# Helpful Sources
- https://guides.github.com/introduction/flow/
- https://www.atlassian.com/git/tutorials
