Param
(
  $InstallPath,
  $ToolsPath,
  $Package,
  $Project
)

# Save current project state.
$Project.Save()

# Load project XML.
$ProjectFullName = $Project.FullName

$Doc = New-Object System.Xml.XmlDocument
$Doc.Load($ProjectFullName)

$Namespace = 'http://schemas.microsoft.com/developer/msbuild/2003'

# Find the node containing the file. The tag "Content" may be replace by "None" depending of the case, check your .csproj file.
$ContentNodes = Select-Xml "//msb:Project/msb:ItemGroup/msb:Content[starts-with(@Include, 'rubberband-dll-') and (substring(@Include, string-length(@Include) - 3) = '.dll')]" $Doc -Namespace @{msb = $Namespace}

#check if the node exists.
If ($ContentNodes -ne $Null)
{
	$CopyNodeName = "CopyToOutputDirectory"

	ForEach ($ContentNode In $ContentNodes)
	{
		$NoneNode = $Doc.CreateElement("None", $Namespace)
		$NoneNode.SetAttribute("Include", $ContentNode.Node.Attributes["Include"].Value)
		
		$CopyNode = $Doc.CreateElement($CopyNodeName, $Namespace)
		$CopyNode.InnerText = "PreserveNewest"
		
		$NoneNode.AppendChild($CopyNode)

		$ContentNode.Node.ParentNode.ReplaceChild($NoneNode, $ContentNode.Node)
	}

	$DTE = $Project.DTE

	$Project.Select("vsUISelectionTypeSelect")

	$DTE.ExecuteCommand("Project.UnloadProject")
	$Doc.Save($ProjectFullName)
	$DTE.ExecuteCommand("Project.ReloadProject")
}
