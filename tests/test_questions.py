# -*- coding: utf-8 -*-
"""
Testing user input stuff is more difficult but this might help:
http://blog.broekhuizen.nu/pytest-mock-builtins.html
"""

from .context import questions


def test_build_prompt():
  # Test simple prompt (100% automatic)
  prompt = questions.build_prompt('name', '(Heisenberg)')
  assert prompt == 'name: (Heisenberg) '

  # Test simple user defined prompt
  prompt = questions.build_prompt("Say my name: {} ", '(Heisenberg)')
  assert prompt == "Say my name: (Heisenberg) "

  # Test a '?' prompt (should not add trailing ':')
  prompt = questions.build_prompt('What is your name?', '(Heisenberg)')
  assert prompt == 'What is your name? (Heisenberg) '

  # Test a more complex 'inline' prompt
  prompt = questions.build_prompt('I {}ed to China', '[walk]')
  assert prompt == 'I [walk]ed to China '
