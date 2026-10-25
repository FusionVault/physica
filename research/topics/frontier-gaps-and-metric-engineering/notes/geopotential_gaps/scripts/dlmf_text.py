"""Render a saved DLMF HTML page as plain text, with each <math> replaced by its LaTeX alttext.
Usage: python3 -I dlmf_text.py page.html > page.txt"""
import sys
from html.parser import HTMLParser

class P(HTMLParser):
    def __init__(self):
        super().__init__()
        self.out, self.math, self.skip = [], 0, 0
    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if tag == 'math':
            if self.math == 0:
                self.out.append(' $' + a.get('alttext', '?') + '$ ')
            self.math += 1
        elif tag in ('script', 'style'):
            self.skip += 1
        elif tag in ('p', 'div', 'tr', 'br', 'h1', 'h2', 'h3', 'h4', 'li', 'table'):
            self.out.append('\n')
    def handle_endtag(self, tag):
        if tag == 'math':
            self.math -= 1
        elif tag in ('script', 'style'):
            self.skip -= 1
    def handle_data(self, d):
        if self.math == 0 and self.skip == 0:
            self.out.append(d)

p = P()
p.feed(open(sys.argv[1], encoding='utf-8', errors='replace').read())
import re
text = ''.join(p.out)
text = re.sub(r'[ \t]+', ' ', text)
text = re.sub(r'\n\s*\n+', '\n', text)
print(text)
