from pyworkflow.gui.browser import FileHandler
import pyworkflow.utils as pwutils


class ImodHandler(FileHandler):

    def getFileActions(self, objFile):
        from imod.viewers import ImodObjectView
        fn = objFile.getPath()
        return [('Open with Imod', lambda: ImodObjectView(fn).show(),
                 pwutils.Icon.ACTION_VISUALIZE)]

    def getFileIcon(self, objFile):
        return pwutils.Icon.FILE_STACK if not objFile.isLink() else pwutils.Icon.FILE_STACK_LINK
