from pyworkflow.gui import FileHandler, pwutils


class ImodHandler(FileHandler):

    def getFileActions(self, objFile):
        from .viewers import ImodObjectView
        fn = objFile.getPath()
        return [('Open with Imod', lambda: ImodObjectView(fn).show(),
                 pwutils.Icon.ACTION_VISUALIZE)]

    def getFileIcon(self, objFile):
        return pwutils.Icon.FILE_STACK if not objFile.isLink() else pwutils.Icon.FILE_STACK_LINK
