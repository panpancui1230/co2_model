import ast

class FunctionAndClassVisitor(ast.NodeVisitor):
    def __init__(self):
        self.functions = []
        self.classes = []

    def visit_FunctionDef(self, node):
        self.functions.append(node.name)
        self.generic_visit(node)

    def visit_ClassDef(self, node):
        class_info = {
            'name': node.name,
            'methods': []
        }
        for elem in node.body:
            if isinstance(elem, ast.FunctionDef):
                class_info['methods'].append(elem.name)
        self.classes.append(class_info)
        self.generic_visit(node)

def extract_functions_and_classes(script):
    tree = ast.parse(script)
    visitor = FunctionAndClassVisitor()
    visitor.visit(tree)
    return visitor.functions, visitor.classes

with open('./KEA3_Q_delete.py', 'r', encoding='utf-8') as file:
    script_content = file.read()

functions, classes = extract_functions_and_classes(script_content)

print(functions, classes)
