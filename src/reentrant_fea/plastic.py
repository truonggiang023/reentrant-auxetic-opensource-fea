"""Read the original plastic table without importing/executing Abaqus."""
import ast
import re
from pathlib import Path

def from_reference(path):
    s=Path(path).read_text(encoding='utf-8-sig')
    m=re.search(r'\bPLASTIC_DATA\s*=\s*\(',s)
    if not m: raise ValueError('PLASTIC_DATA table not found')
    start=s.index('(',m.start())
    depth=0
    end=None
    for i in range(start,len(s)):
        if s[i]=='(': depth+=1
        elif s[i]==')':
            depth-=1
            if depth==0:
                end=i+1
                break
    if end is None: raise ValueError('Unclosed PLASTIC_DATA')
    values=ast.literal_eval(s[start:end])
    if any(len(p)!=2 for p in values): raise ValueError('Invalid plastic table')
    if any(values[i+1][1]<=values[i][1] for i in range(len(values)-1)):
        raise ValueError('Plastic strains must increase')
    return values
