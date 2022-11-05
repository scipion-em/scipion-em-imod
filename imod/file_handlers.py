from pyworkflow.gui import FileHandler, pwutils


class ImodHandler(FileHandler):

    def getFileActions(self, objFile):
        from .viewers import ImodObjectView
        fn = objFile.getPath()
        return [('Open with Imod', lambda: ImodObjectView(fn).show(),
                 pwutils.Icon.ACTION_VISUALIZE)]

    def getFileIcon(self, objFile):
        return 'file_stack.gif' if not objFile.isLink() else 'file_stack_link.gif'
