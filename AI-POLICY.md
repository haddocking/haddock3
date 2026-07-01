# AI Policy

Using AI tools to help write code is fine.

That said, HADDOCK3 is a scientific software — a mistake in a scoring function or restraint definition can silently corrupt research results.

## What's fine

Using Copilot, ChatGPT, Claude or any similar AI agent to draft tests, improve documentation, write boilerplate, or understand code before modifying it is all fine.

We care about the quality of what ends up in the codebase, not how you got there.

## What isn't fine

Don't trust AI output for anything that requires scientific judgment: scoring functions, energy terms, restraint logic, docking algorithms, CNS/topology parameter handling. 

If AI substantially shaped any of that kind of work, disclose this in the PR and explain how you verified it.

## You own what you submit

- Read the diff.
- Understand every line.
- Be ready to explain and/or discuss proposed changes in review.

A reviewer might ask questions about your code, if you can't answer it because the AI wrote that part and you do not understand why, the PR might be closed.

PR reviews, descriptions and responses should be written by yourself.

## Disclosure

If AI did a meaningful chunk of the work, mention it briefly in the PR description — something like "drafted with Copilot, verified manually" is enough. You don't need to declare autocomplete or grammar fixes.

--- 

If you have questions or concerns, open an issue or send an e-mail to `bonvinlab.support@uu.nl`.
