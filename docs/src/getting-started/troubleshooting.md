---
title: Troubleshooting - AUSAXS
description: Solutions to common problems when installing or running AUSAXS on Windows, macOS and Linux.
---

# Troubleshooting

Getting the programs to run on some platforms can be somewhat tricky. Please refer to the following if you encounter any issues. 

> [!NOTE]
> Can't find your issue? Please create a [discussion thread](https://github.com/AUSAXS/AUSAXS/discussions) with more information so we can help get it solved. 

## Windows

<details markdown="1">
<summary>Why do I get SmartScreen warnings?</summary>
The Microsoft SmartScreen filter may warn you about potential security risks when running AUSAXS for the first time. You should be able to ignore this warning by clicking "More information", where you will be able to override the warning. Removing this notice from our end requires a yearly fee to sign our program, which we are currently unable to support. <br><br>
  
**Warning!** You should only ignore SmartScreen warnings for software you trust! AUSAXS is completely open source, so if you have any security concerns, you are free to check the code. 
</details>


## MacOS

<details markdown="1">
<summary>Why does it immediately crash upon opening the program?</summary>

The binaries are compiled against specific MacOS versions, which are usually not backwards compatible. We strive to compile against the [oldest MacOS versions available](https://github.com/actions/runner-images?tab=readme-ov-file#available-images) on GitHub. <br><br>

In case they do not work on your system, you can try to compile the software yourself as described [here](compiling.md#macos).
</details>

## Linux
Usually there are no issues on this platform. Perhaps you simply lack the permissions to execute the program? Try marking them as executable: 
`chmod +x <program>`
