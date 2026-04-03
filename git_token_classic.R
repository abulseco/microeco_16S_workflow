# Get classic token from github
# See here for more info: https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/managing-your-personal-access-tokens

# Install the `usethis` and `gitcreds` packages
install.packages(c("usethis", "gitcreds"))

# Add your GitHub username and email
usethis::use_git_config(user.name = "XX",
                        user.email = "XX") # replace with your own Github info

# Create a token (Note this will open a GitHub browser tab)
## See steps 6-10 in GitHub's PAT tutorial (link below)
usethis::create_github_token()

# Now, give your token to RStudio
## After you run this line follow the prompts in the "Console" tab of RStudio
gitcreds::gitcreds_set()