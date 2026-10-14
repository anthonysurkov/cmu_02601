'''
Node Class (can be changed if needed):

Attributes:
- sequence (str)
- children (list[Node | str])
'''


from __future__ import annotations


class Node:
    def __init__(self, sequence: str, children: list[Node | str] | None = None):
        self.sequence = sequence
        self.children = children if children is not None else []
