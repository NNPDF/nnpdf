from validphys.app import App
from validphys import comparefittemplates


class CompareFitApp(App):
    def add_positional_arguments(self, parser):
        pass

    def get_commandline_arguments(self, cmdline=None):
        args = super().get_commandline_arguments(cmdline)
        runcard = comparefittemplates.template_path
        # This is needed because the environment wants to know how to resolve
        # the relative paths to find the templates. Best to have the template
        # look as much as possible as a runcard passed from the command line
        args['config_yml'] = runcard
        return args


def main():
    a = CompareFitApp()
    a.main()

if __name__ == '__main__':
    main()
