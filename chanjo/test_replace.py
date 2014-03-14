#!/usr/bin/env python

from clint.textui import puts, colored


def input_replace(prompt, default=None, color=colored.green):
  """Use in place in input. Asks a question, waits for user input, updates the
  same line by replacing the default.

  .. code-block::

    >>> my_name = input_replace('Say my name: {} ', 'Heisenberg')
    # Wait for user input
    Say my name: (Heisenberg) Walter
    # Updates **the same** line with 'Walter' in green
    Say my name: Walter
    >>> print(my_name)
    Walter

  Inspired by Bower's init prompt which confirms user input by replacing the
  default option in line. Ref: http://stackoverflow.com/questions/12586601

  Args:
    prompt (str): Question to print to user, '{}' will be replaced by default
    default (str, optional): Default option unless replaced by user
    color (function, optional): :func:`clint.textui.colored` function

  Returns:
    str: User input or default
  """
  # Helper variables
  CURSOR_UP_ONE = '\x1b[1A'
  ERASE_LINE = '\x1b[2K'

  # Determine if a default was submitted
  if default:
    # Prepare the default-part of the prompt
    default_string = '({})'.format(default)
  else:
    # Not relevant since ``promt`` shouldn't include a '{}'
    default_string = ''

  # Pass question to user and wait for response
  # Write default option in parentheses, use it as response if nothing
  # was submitted by user.
  response = input(prompt.format(default_string)) or default

  # Print the updated confirmation line by replacing the previous
  print(CURSOR_UP_ONE + ERASE_LINE + prompt.format(color(response)))

  return response

my_name = input_replace('Say my name: {} ', 'Heisenberg')
