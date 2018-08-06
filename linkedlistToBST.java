/**
 * Definition for singly-linked list.
 * public class ListNode {
 *     int val;
 *     ListNode next;
 *     ListNode(int x) { val = x; }
 * }
 */
/**
 * Definition for a binary tree node.
 * public class TreeNode {
 *     int val;
 *     TreeNode left;
 *     TreeNode right;
 *     TreeNode(int x) { val = x; }
 * }
 */
class Solution {
    public TreeNode sortedListToBST(ListNode head) {
        if(head==null)
            return null;
        
        List<Integer> valList=new ArrayList<Integer>();
        while(head!=null){
            valList.add(head.val);
            head=head.next;
        }
        return binaryHelper(0,valList.size()-1,valList);
    }
    private TreeNode binaryHelper(int strt, int end, List<Integer> lst){
        if(strt>end)
            return null;
        if(strt==end){
            return new TreeNode(lst.get(strt));
        }
        
        
        int split=(strt+end)/2;
        TreeNode root=new TreeNode(lst.get(split));
        root.left=binaryHelper(strt,split-1,lst);
        root.right=binaryHelper(split+1,end,lst);
        return root;
    }
}
