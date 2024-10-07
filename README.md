
TEM-Collagen-Analysis-Plugin

This is a custom plugin developed for ImageJ version 1. The plugin allows users to analyze clusters of collagen fibrils. This guide will help you clone the repository, build the plugin, and install it into your local ImageJ installation.

Prerequisites
Before proceeding, make sure you have the following installed:

ImageJ version 1 - Download and install from ImageJ website.
Java Development Kit (JDK) - Version 8 or higher.
Maven - A Java build tool, you can install it from Maven's official site.
Instructions to Clone, Build, and Install
1. Clone the Repository
First, clone the plugin repository to your local machine. Open your terminal (or command prompt) and run:

https://github.com/Ty-Watson/TEM-Collagen-Analysis-Plugin.git
This will download all the files from the repository into a local directory.

2. Navigate to the Plugin Directory
Change to the directory where the plugin's source code is located:

cd your-plugin-repo

3. Build the JAR File
Use Maven to build the plugin into a JAR file. This will compile the Java files and package them into a plugin JAR file:

mvn clean package

If everything is set up correctly, this command will create a .jar file in the target directory.

4. Install the Plugin
Option 1: Manual Installation
Once the JAR file is created, you need to copy it into ImageJ’s plugins folder:

Locate your ImageJ installation directory.
Navigate to the plugins folder inside your ImageJ directory.
On Windows, this could be something like C:\Program Files\ImageJ\plugins.
On macOS/Linux, it could be /Applications/ImageJ/plugins or ~/ImageJ/plugins.
Copy the JAR file (your-plugin-repo.jar) from the target folder into the plugins folder.
Option 2: Use Maven for Installation
Die-hard command-line developers can use Maven directly by calling mvn in the project root.

However you build the project, in the end, you will have the .jar file (called an artifact in Maven speak) in the target/ subdirectory.

To copy the artifact into the correct place, you can call:


mvn -Dscijava.app.directory=/path/to/ImageJ.app/

This will not only copy your artifact but also all the dependencies. Restart your ImageJ or call Help › Refresh Menus to see your plugin in the menus.

5. Run ImageJ and Use the Plugin
Open ImageJ.
You should see your plugin listed under the Plugins menu in ImageJ.
You can now use the plugin by selecting it from the Plugins menu!

6. Updating the Plugin
To update the plugin, pull the latest changes from the repository and rebuild the JAR file:

git pull
mvn clean package

Then, copy the new JAR file to your ImageJ plugins directory using either the manual method or the Maven method and restart ImageJ.

Troubleshooting
Build Errors: Ensure you have the correct version of Java and Maven installed, and that your environment variables (JAVA_HOME and MAVEN_HOME) are properly set.
Plugin Not Showing: Ensure the JAR file is placed in the correct plugins folder or installed via Maven and restart ImageJ.
If you encounter any issues, feel free to open an issue on the GitHub repository.

License
This project is licensed under the MIT License. See the LICENSE file for details.

Author
[Your Name] - Developer of the plugin.

This addition will guide users on how to use Maven to directly copy the plugin and its dependencies into the ImageJ directory, as requested.






