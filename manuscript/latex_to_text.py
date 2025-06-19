#!/usr/bin/env python3
"""
Simple LaTeX to text converter for previewing manuscript content
"""

import re
import sys

def latex_to_text(latex_content):
    """Convert LaTeX content to readable text"""
    
    # Remove comments
    text = re.sub(r'%.*$', '', latex_content, flags=re.MULTILINE)
    
    # Remove document class and packages
    text = re.sub(r'\\documentclass.*?\n', '', text)
    text = re.sub(r'\\usepackage.*?\n', '', text)
    text = re.sub(r'\\geometry.*?\n', '', text)
    
    # Remove preamble commands
    text = re.sub(r'\\title\{([^}]*)\}', r'TITLE: \1', text)
    text = re.sub(r'\\author\{([^}]*)\}', r'AUTHORS: \1', text)
    text = re.sub(r'\\date\{([^}]*)\}', r'DATE: \1', text)
    
    # Convert sections
    text = re.sub(r'\\section\{([^}]*)\}', r'\n\n=== \1 ===\n', text)
    text = re.sub(r'\\subsection\{([^}]*)\}', r'\n\n--- \1 ---\n', text)
    text = re.sub(r'\\subsubsection\{([^}]*)\}', r'\n\n\1:\n', text)
    
    # Remove document structure
    text = re.sub(r'\\begin\{document\}', '', text)
    text = re.sub(r'\\end\{document\}', '', text)
    text = re.sub(r'\\maketitle', '', text)
    
    # Handle abstracts
    text = re.sub(r'\\begin\{abstract\}', '\n=== ABSTRACT ===\n', text)
    text = re.sub(r'\\end\{abstract\}', '\n', text)
    
    # Remove figures and tables (keep captions)
    text = re.sub(r'\\begin\{figure\}.*?\\caption\{([^}]*)\}.*?\\end\{figure\}', 
                  r'\n[FIGURE: \1]\n', text, flags=re.DOTALL)
    text = re.sub(r'\\begin\{table\}.*?\\caption\{([^}]*)\}.*?\\end\{table\}', 
                  r'\n[TABLE: \1]\n', text, flags=re.DOTALL)
    
    # Remove other environments we can't easily convert
    text = re.sub(r'\\begin\{[^}]+\}', '', text)
    text = re.sub(r'\\end\{[^}]+\}', '', text)
    
    # Handle citations
    text = re.sub(r'\\cite\{([^}]*)\}', r'[REF: \1]', text)
    text = re.sub(r'\\citep\{([^}]*)\}', r'(\1)', text)
    text = re.sub(r'\\citet\{([^}]*)\}', r'\1', text)
    
    # Handle emphasis
    text = re.sub(r'\\textbf\{([^}]*)\}', r'**\1**', text)
    text = re.sub(r'\\textit\{([^}]*)\}', r'*\1*', text)
    text = re.sub(r'\\emph\{([^}]*)\}', r'*\1*', text)
    
    # Handle special characters
    text = re.sub(r'\\\&', '&', text)
    text = re.sub(r'\\%', '%', text)
    text = re.sub(r'\\\$', '$', text)
    
    # Remove remaining LaTeX commands
    text = re.sub(r'\\[a-zA-Z]+(\[[^\]]*\])?\{[^}]*\}', '', text)
    text = re.sub(r'\\[a-zA-Z]+', '', text)
    
    # Clean up whitespace
    text = re.sub(r'\n\s*\n\s*\n', '\n\n', text)
    text = re.sub(r'^\s+', '', text, flags=re.MULTILINE)
    
    return text.strip()

def main():
    """Main function to convert LaTeX file to text"""
    
    # Read the LaTeX file
    try:
        with open('soybean_drought_rnaseq.tex', 'r', encoding='utf-8') as f:
            latex_content = f.read()
    except FileNotFoundError:
        print("Error: soybean_drought_rnaseq.tex not found!")
        sys.exit(1)
    
    # Convert to text
    text_content = latex_to_text(latex_content)
    
    # Write to output file
    output_file = 'manuscript_preview.txt'
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(text_content)
    
    print(f"LaTeX content converted to readable text: {output_file}")
    print(f"Text length: {len(text_content)} characters")
    
    # Show first few lines
    print("\n=== PREVIEW (first 20 lines) ===")
    lines = text_content.split('\n')
    for i, line in enumerate(lines[:20]):
        print(f"{i+1:2d}: {line}")
    
    if len(lines) > 20:
        print(f"... and {len(lines) - 20} more lines")

if __name__ == "__main__":
    main() 