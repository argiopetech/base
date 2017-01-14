########################################
II. Skip the install and go to the cloud
########################################

BASE-9 executables are available as Amazon Elastic Compute Cloud (EC2) images. Up-to-date instance IDs are listed in the release descriptions at http://github.com/argiopetech/base/releases. An Amazon Web Services (AWS) account is required to use these instances.

To run your code on EC2,

1. Log in to your AWS account
2. Navigate to EC2
3. Navigate to the “Instances” pane
4. Click “Launch Instance”
5. Choose “Community AMIs”
6. Enter the AMI code for the version of BASE-9 you would like to run
7. Select “Review and Launch”, then “Launch”
8. Wait for the instance to launch
9. use rsync/scp to copy your data to the public IP of the instance
10. Login to your instance with SSH

The default user name is ``ec2-user``. There is no root password by default.

BASE-9 executables are in ``/usr/local/bin`` (should be in the path). Current models
(appropriate for the installed version of BASE-9) are in ``/usr/local/share/base-models``.

The instance operating system is the newest release of FreeBSD 10. The tcsh, csh, sh, and bash
shells are available. VI, VIM, Emacs, and nano are pre-installed.